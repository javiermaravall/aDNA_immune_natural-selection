#!/usr/bin/env python3
import os
import argparse
import numpy as np
import polars as pl
from collections import defaultdict

# -------------------- CONSTANT PATHS (v5) --------------------
BASELINELD_FMT = "ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD.{chr}.annot.gz"
CV2F_FMT       = "cV2F_annot/cV2F.{chr}.annot.gz"

RSID_CORR_FMT  = "hg19-ID_to_rsID_correspondence_chr{chr}.tsv"
BLOCKS_FMT     = "jackknife_200_blocks/blocks_chr{chr}.tsv"

LMM_TABLE      = "LMM_withZ_GLMM.corrected.tsv"
FM_DIR         = "fine-mapping/Selection/results"
FM_FILE_FMT    = FM_DIR + "/fine-mapping_results_final_chr{chr}.tsv"

OUT_BASE = "functional_enrichment/raw_results"
os.makedirs(OUT_BASE, exist_ok=True)
os.makedirs(os.path.join(OUT_BASE, "results_by_block"), exist_ok=True)

# Three PIP ranges (keys used in column names)
RANGE_KEYS = [">0.95", "0.5-0.95", "0.05-0.5"]

# -------------------- ARGPARSE --------------------
def parse_args():
    p = argparse.ArgumentParser(description="BaselineLD/cV2F functional enrichment with jackknife (v5, 3 ranges).")
    p.add_argument("category", type=str, help="Annotation column: e.g., coding_UCSC, conservedLindbladToh, or cV2F_binarized")
    return p.parse_args()

args = parse_args()
category = args.category
print(f"[INFO] Processing annotation category: {category}")

# -------------------- HELPERS --------------------
def read_baselineld_or_cv2f(chrom: int, category: str) -> pl.DataFrame:
    """
    Return DF with columns: rsID, <category>.
    If category == 'cV2F_binarized', load from CV2F_FMT and pick that column.
    Otherwise, load from BASELINELD_FMT and pick the requested baselineLD column.
    """
    if category == "cV2F_binarized":
        file_path = CV2F_FMT.format(chr=chrom)
        col_needed = "cV2F_binarized"
    else:
        file_path = BASELINELD_FMT.format(chr=chrom)
        col_needed = category

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Annotation file not found: {file_path}")

    df = pl.read_csv(file_path, separator="\t", infer_schema_length=100000)
    if "SNP" not in df.columns:
        raise ValueError(f"'SNP' column missing in {file_path}")
    if col_needed not in df.columns:
        raise ValueError(f"Column '{col_needed}' not found in {file_path}")

    return df.select([pl.col("SNP").alias("rsID"), pl.col(col_needed).alias(category)])

def read_correspondence_table(chrom: int) -> pl.DataFrame:
    """Return DF with columns: ID (Ali hg19), rsID."""
    path = RSID_CORR_FMT.format(chr=chrom)
    if not os.path.exists(path):
        raise FileNotFoundError(f"rsID correspondence file not found: {path}")
    df = pl.read_csv(path, separator="\t").rename({"ID_Ali_hg19": "ID"})
    if not {"ID", "rsID"}.issubset(df.columns):
        raise ValueError(f"Expected 'ID_Ali_hg19' and 'rsID' in {path}")
    return df.select(["ID", "rsID"])

def read_blocks(chrom: int) -> dict[str, list[str]]:
    """Return dict block_id -> list of Ali IDs (variant_id)."""
    block_file = BLOCKS_FMT.format(chr=chrom)
    if not os.path.exists(block_file):
        raise FileNotFoundError(f"Block file not found: {block_file}")
    df = pl.read_csv(block_file, separator="\t")
    if not {"block_id", "variant_id"}.issubset(df.columns):
        raise ValueError(f"Block file missing required columns in {block_file}")
    blocks = defaultdict(list)
    for block_id, variant_id in df.select(["block_id", "variant_id"]).iter_rows():
        blocks[block_id].append(variant_id)
    return blocks

def assign_category_join(df_ids: pl.DataFrame, corr: pl.DataFrame, annot: pl.DataFrame, chrom: int, category: str, label: str) -> pl.DataFrame:
    """
    Vectorized assignment:
      ID --(join)--> rsID --(join)--> annotation value  => Category bool or null.
    Returns DF with columns: ID, (optionally PIP), Category (bool).
    """
    df = (
        df_ids.join(corr, on="ID", how="left")
              .join(annot, on="rsID", how="left")
              .with_columns(
                  pl.when(pl.col(category).is_null())
                    .then(pl.lit(None))
                    .otherwise(pl.col(category) != 0)
                    .alias("Category")
              )
              .filter(pl.col("Category").is_not_null())
    )
    keep = [c for c in ["ID", "PIP", "Category"] if c in df.columns]
    out = df.select(keep)
    print(f"  [INFO] chr{chrom} {label}: total={len(df_ids)}; assigned={len(out)}; unassigned={len(df_ids) - len(out)}")
    return out

# -------------------- LOAD LMM ID UNIVERSE ONCE --------------------
lmm_all = (
    pl.scan_csv(LMM_TABLE, separator="\t")
      .select(["CHROM", "ID"])   # Ali hg19 ID
      .collect()
)
print(f"[INFO] Loaded LMM ID universe: n={len(lmm_all)}")

# -------------------- MAIN --------------------
rsid_data_parts = []
variant_data_parts = []
all_blocks = defaultdict(list)

for num in range(1, 23):
    print(f"\n[INFO] Processing chromosome {num}")

    # Blocks (IDs)
    blocks_chr = read_blocks(num)
    for block_id, ids in blocks_chr.items():
        all_blocks[block_id].extend(ids)

    # Preload correspondence & annotation once per chromosome
    corr = read_correspondence_table(num)              # ID, rsID
    annot = read_baselineld_or_cv2f(num, category)     # rsID, <category>

    # LMM IDs for this chromosome (universe for p2_category)
    lmm_chr_ids = lmm_all.filter(pl.col("CHROM") == num).select(["ID"])
    rsid_df_chr = assign_category_join(lmm_chr_ids, corr, annot, num, category, label="LMM")
    rsid_data_parts.append(rsid_df_chr.select(["ID", "Category"]))

    # Fine-mapping for this chromosome
    fm_path = FM_FILE_FMT.format(chr=num)
    if os.path.exists(fm_path):
        fm = (
            pl.read_csv(fm_path, separator="\t")
              .select(["ID", "PIP_coloc"])
              .rename({"PIP_coloc": "PIP"})
        )
        fm_with_cat = assign_category_join(fm, corr, annot, num, category, label="FM")
        variant_data_parts.append(fm_with_cat)
    else:
        print(f"  [WARN] Missing fine-mapping file: {fm_path}")

# Concatenate all chromosomes
rsid_data = pl.concat(rsid_data_parts) if rsid_data_parts else pl.DataFrame({"ID": [], "Category": []})
variant_data = pl.concat(variant_data_parts) if variant_data_parts else pl.DataFrame({"ID": [], "PIP": [], "Category": []})

print(f"\n[INFO] rsid_data (assigned) n={len(rsid_data)}")
print(f"[INFO] variant_data (assigned) n={len(variant_data)}\n")

# -------------------- p2_category (fractions over assigned rsIDs) -----------
# Category is boolean (True/False)
cat_counts = rsid_data.group_by("Category").len()
total_assigned = int(len(rsid_data))
p2_category = {row["Category"]: (row["len"] / total_assigned if total_assigned > 0 else 0.0)
               for row in cat_counts.iter_rows(named=True)}
print(f"[INFO] p2_category = {p2_category}\n")

# -------------------- Build ID sets for 3-range jackknife speed -------------
# rsid universe per Category
rsid_ids_by_cat = {
    cv: set(rsid_data.filter(pl.col("Category") == cv).select("ID").to_series().to_list())
    for cv in [True, False]
}

# Three PIP ranges (exact boundaries)
range_exprs = {
    ">0.95":    (pl.col("PIP") > 0.95),
    "0.5-0.95": (pl.col("PIP") >= 0.5) & (pl.col("PIP") <= 0.95),
    "0.05-0.5": (pl.col("PIP") >= 0.05) & (pl.col("PIP") < 0.5),
}

# For each range, collect unique IDs in category True/False
range_sets_by_cat = {True: {}, False: {}}
for rk, expr in range_exprs.items():
    df_r = variant_data.filter(expr).select(["ID", "Category"])
    for cv in [True, False]:
        ids = df_r.filter(pl.col("Category") == cv).select("ID").to_series().to_list()
        range_sets_by_cat[cv][rk] = set(ids)

# For enrichment, per-range denominator is the size of the PIP subset = sum over categories
fr_totals = {rk: sum(len(range_sets_by_cat[cv][rk]) for cv in [True, False]) for rk in RANGE_KEYS}

# Fraction of variants in PIP subset that belong to each category
fr_results = {
    cv: {rk: (len(range_sets_by_cat[cv][rk]) / fr_totals[rk] if fr_totals[rk] > 0 else 0.0)
         for rk in RANGE_KEYS}
    for cv in [True, False]
}

# Baseline enrichment
enrichment_baseline = {
    cv: {f"Enrichment_{rk}": (fr_results[cv][rk] / p2_category.get(cv, 0.0) if p2_category.get(cv, 0.0) > 0 else 0.0)
         for rk in RANGE_KEYS}
    for cv in [True, False]
}

print(f"[INFO] fr_totals = {fr_totals}")
print(f"[INFO] fr_results = {fr_results}")
print(f"[INFO] enrichment (baseline) = {enrichment_baseline}\n")

# -------------------- Jackknife (set-based, 3 ranges) -----------------------
jackknife_results = defaultdict(lambda: defaultdict(list))
results_by_block = defaultdict(list)

total_blocks = len(all_blocks)
progress_file = os.path.join(OUT_BASE, f"progress_{category}.txt")

print(f"[INFO] Jackknife across {total_blocks} blocks\n")
with open(progress_file, "w") as progress:
    for i, (block_id, block_variant_ids) in enumerate(all_blocks.items(), 1):
        block_set = set(block_variant_ids)

        # p2_category excluding block (over assigned universe)
        jk_counts_by_cat = {cv: len(rsid_ids_by_cat[cv] - block_set) for cv in [True, False]}
        jk_total = sum(jk_counts_by_cat.values())
        jk_p2_category = {cv: (jk_counts_by_cat[cv] / jk_total) if jk_total > 0 else 0.0 for cv in [True, False]}

        # PIP-range totals excluding block = sum over categories
        jk_fr_totals = {
            rk: sum(len((range_sets_by_cat[cv][rk] - block_set)) for cv in [True, False])
            for rk in RANGE_KEYS
        }

        # Fractions within each PIP subset
        jk_fr_results = {
            cv: {rk: (len(range_sets_by_cat[cv][rk] - block_set) / jk_fr_totals[rk] if jk_fr_totals[rk] > 0 else 0.0)
                 for rk in RANGE_KEYS}
            for cv in [True, False]
        }

        # Enrichment excluding block
        jk_enrichment = {
            cv: {f"Enrichment_{rk}": (jk_fr_results[cv][rk] / jk_p2_category[cv] if jk_p2_category[cv] > 0 else 0.0)
                 for rk in RANGE_KEYS}
            for cv in [True, False]
        }

        # Accumulate results
        for cv in [True, False]:
            for rk in RANGE_KEYS:
                jackknife_results[cv][f"Enrichment_{rk}"].append(jk_enrichment[cv][f"Enrichment_{rk}"])

        # Per-block rows (one row per category per block)
        for cv in [True, False]:
            results_by_block[block_id].append({
                "Category": cv,
                "Enrichment_>0.95":   jk_enrichment[cv]["Enrichment_>0.95"],
                "Enrichment_0.5-0.95": jk_enrichment[cv]["Enrichment_0.5-0.95"],
                "Enrichment_0.05-0.5": jk_enrichment[cv]["Enrichment_0.05-0.5"],
            })

        # Progress
        progress_message = f"{i} of {total_blocks} blocks processed — {i / total_blocks * 100:.2f}% done\n"
        print(progress_message.strip())
        progress.write(progress_message)
        progress.flush()

print("\n[INFO] Jackknife complete\n")

# -------------------- Aggregate jackknife (mean & SE) -----------------------
final_results = {cv: {} for cv in [True, False]}
for cv in [True, False]:
    for rk in RANGE_KEYS:
        arr = jackknife_results[cv][f"Enrichment_{rk}"]
        mean_est = float(np.mean(arr)) if arr else 0.0
        # Match your earlier SE formula:
        se_est   = float(np.std(arr) * np.sqrt(len(arr) - 1)) if len(arr) > 1 else 0.0
        final_results[cv][f"Enrichment_{rk}"]    = mean_est
        final_results[cv][f"Enrichment_{rk}_SE"] = se_est

# -------------------- Write summary (3 ranges) ------------------------------
out_file = os.path.join(OUT_BASE, f"enrichment_jackknife_by-range_200-blocks_{category}.tsv")
with open(out_file, "w") as f:
    f.write("category\tfr_overall\t"
            "fr_>0.95\tfr_0.5-0.95\tfr_0.05-0.5\t"
            "Enrichment_>0.95\tEnrichment_>0.95_SE\t"
            "Enrichment_0.5-0.95\tEnrichment_0.5-0.95_SE\t"
            "Enrichment_0.05-0.5\tEnrichment_0.05-0.5_SE\n")
    for cv in [True, False]:
        f.write(
            f"{cv}\t{p2_category.get(cv, 0.0)}\t"
            f"{fr_results[cv]['>0.95']}\t{fr_results[cv]['0.5-0.95']}\t{fr_results[cv]['0.05-0.5']}\t"
            f"{final_results[cv]['Enrichment_>0.95']}\t{final_results[cv]['Enrichment_>0.95_SE']}\t"
            f"{final_results[cv]['Enrichment_0.5-0.95']}\t{final_results[cv]['Enrichment_0.5-0.95_SE']}\t"
            f"{final_results[cv]['Enrichment_0.05-0.5']}\t{final_results[cv]['Enrichment_0.05-0.5_SE']}\n"
        )

# -------------------- Write per-block file (3 ranges) -----------------------
by_block_file = os.path.join(OUT_BASE, "results_by_block", f"results_by_block_{category}.tsv")
with open(by_block_file, "w") as f:
    f.write("block_id\tcategory\t"
            "Enrichment_>0.95\tEnrichment_0.5-0.95\tEnrichment_0.05-0.5\n")
    for block_id, rows in results_by_block.items():
        for r in rows:
            f.write(f"{block_id}\t{r['Category']}\t"
                    f"{r['Enrichment_>0.95']}\t"
                    f"{r['Enrichment_0.5-0.95']}\t"
                    f"{r['Enrichment_0.05-0.5']}\n")

print(f"[OK] Wrote summary: {out_file}")
print(f"[OK] Wrote per-block: {by_block_file}")
print("[INFO] Analysis completed.")
