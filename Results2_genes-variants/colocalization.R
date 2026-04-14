#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# ------------------ Thresholds ------------------
p_value_threshold <- 5e-7
#p_value_threshold <- 1e-5
chi_sq_threshold  <- qchisq(1 - p_value_threshold, df = 1)

# ------------------ Helpers (aligned with GWAS script) ------------------
colname_get <- function(dt, target) {
  nm <- names(dt)
  idx <- match(tolower(target), tolower(nm))
  if (!is.na(idx)) nm[idx] else NULL
}
# Try candidates (case-insensitive); first present is used
colname_get_any <- function(dt, candidates) {
  for (cand in candidates) {
    hit <- colname_get(dt, cand)
    if (!is.null(hit)) return(hit)
  }
  NULL
}
all_prop_missing <- function(vec) {
  if (is.null(vec)) return(TRUE)
  vchar <- as.character(vec)
  all(is.na(vchar) | vchar %in% c("NA", "NA_sdY", ""))
}
stopf <- function(...) stop(sprintf(...), call. = FALSE)

# ------------------ Args ------------------
# Usage: coloc_gwas_like.R <trait> [N_for_Trait2_optional]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: coloc_gwas_like.R <trait> [N_for_Trait2_optional]", call. = FALSE)
}
trait2_name <- args[1]
N_trait2_arg <- if (length(args) >= 2) as.numeric(args[2]) else NA_real_
if (!is.na(N_trait2_arg) && (N_trait2_arg <= 0 || is.nan(N_trait2_arg))) stop("N_for_Trait2 must be positive if provided")
cat(sprintf("[INFO] Trait2 (GWAS trait): %s | N_arg=%s\n", trait2_name, ifelse(is.na(N_trait2_arg), "NA", as.character(N_trait2_arg))))

# ------------------ Paths ------------------
outdir <- sprintf(
  "colocalization/GWAS/results_by_trait/%s",
  trait2_name
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

windows_file   <- "fine-mapping/windows/windows.tsv"
selection_file <- "LMM_withZ_GLMM.corrected.tsv"

# GWAS trait path (same convention as your GWAS script)
gwas_file_tpl <- "data/1.GWAS/all_sources/%s.tsv.gz"
gwas_file <- sprintf(gwas_file_tpl, trait2_name)
if (!file.exists(gwas_file)) stop(sprintf("GWAS trait file not found: %s", gwas_file))

# ------------------ Main ------------------

# Read windows once (unchanged)
stopifnot(file.exists(windows_file))
windows_all <- fread(windows_file)

# Loop autosomes
for (num in 1:22) {
  cat(sprintf("\n[CHR %d] Starting...\n", num))

  # Windows for this chromosome
  windows_dt <- copy(windows_all)
  windows_dt[, chr := as.integer(sub("^chr([0-9]+)_.*$", "\\1", window_ID))]
  windows_dt <- windows_dt[chr == num]
  windows_dt[, chr := NULL]
  if (nrow(windows_dt) == 0L) { cat(sprintf("[CHR %d] No windows. Skip.\n", num)); next }

  # --------------------------
  # Read entire GWAS chromosome once
  # --------------------------
  cmd_chr <- sprintf("gunzip -c %s | awk -F '\\t' 'NR==1 || ($3==%d)'", shQuote(gwas_file), num)
  gwas_chr <- tryCatch(fread(cmd = cmd_chr), error = function(e) {
    cat(sprintf("[CHR %d] Failed to read gwas chr with cmd: %s\nError: %s\n", num, cmd_chr, e$message))
    return(NULL)
  })
  if (is.null(gwas_chr) || nrow(gwas_chr) == 0L) { cat(sprintf("[CHR %d] No GWAS chr rows. Skip.\n", num)); next }

  # Normalize column names for GWAS file (aligned with GWAS script)
  col_ID   <- colname_get_any(gwas_chr, c("ID_hg19"))
  col_CHROM<- colname_get_any(gwas_chr, c("CHROM","CHR"))
  col_POS  <- colname_get_any(gwas_chr, c("POS","POSITION"))
  col_BETA <- colname_get_any(gwas_chr, c("BETA","Effect","Beta"))
  col_SE   <- colname_get_any(gwas_chr, c("SE","SE_beta","StdErr","SE_BETA"))
  col_PVAL <- colname_get_any(gwas_chr, c("P_value","PVAL","PVALUE","p","pvalue","P"))
  col_Z    <- colname_get_any(gwas_chr, c("Z","z"))
  col_MAF  <- colname_get_any(gwas_chr, c("MAF","Freq1","AF","EAF","ALT_AF"))
  col_N    <- colname_get_any(gwas_chr, c("N","N_total","Ntot","Ncase+Ncontrol","N_samples"))
  col_Prop <- colname_get_any(gwas_chr, c("Prop_cases","PropCases","prop_cases","prop"))

  # Precompute medN (chromosome level)
  medN_chr <- NA_real_
  if (!is.null(col_N)) {
    Nvals_chr <- suppressWarnings(as.numeric(gwas_chr[[col_N]]))
    medN_chr <- median(Nvals_chr, na.rm = TRUE)
    if (is.nan(medN_chr)) medN_chr <- NA_real_
  }

  # Read selection data for this chromosome (unchanged)
  cmd_sel_chr <- sprintf("awk -F '\\t' 'NR==1 || ($1==%d)' %s", num, shQuote(selection_file))
  selection_chr <- fread(cmd = cmd_sel_chr)
  if ("POS" %in% names(selection_chr)) selection_chr[, POS := as.integer(POS)]
  # Detect selection z column (use z_GLMM.corrected as requested)
  col_Z_sel <- colname_get_any(selection_chr, c("z_GLMM.corrected","z_glmm.corrected","z_GLMM_corrected","z_glmm_corrected"))

  # Collectors
  snp_by_window_list <- list()
  pp_by_window_list  <- list()

  # Iterate windows (selection part remains as in your script)
  for (i in seq_len(nrow(windows_dt))) {
    window <- windows_dt[i]
    start <- window$start
    end <- window$end
    start_central <- window$start_central_1Mb
    end_central <- window$end_central_1Mb
    window_ID <- window$window_ID

    # ---- Selection slice (unchanged slicing) ----
    selection_data <- selection_chr[POS >= start & POS < end]
    if (nrow(selection_data) == 0L) {
      cat(sprintf("\r[CHR %d] Window %d/%d: empty selection (skip)", num, i, nrow(windows_dt)))
      next
    }

    # ---- GWAS Trait2 slice (same approach as GWAS script) ----
    if (is.null(col_POS)) stopf("Missing POS column overall; cannot proceed (chr %d, window %s).", num, window_ID)
    pos_vec_chr <- as.integer(gwas_chr[[col_POS]])
    mask_win <- (pos_vec_chr >= start) & (pos_vec_chr < end)
    gwas_window <- gwas_chr[mask_win, , drop = FALSE]
    if (nrow(gwas_window) == 0L) {
      cat(sprintf("\r[CHR %d] Window %d/%d: empty Trait2 (skip)", num, i, nrow(windows_dt)))
      next
    }

    # ---- CENTRAL significance checks (UPDATED) ----
    # Selection: use z_GLMM.corrected (required)
    if (is.null(col_Z_sel) || !(col_Z_sel %in% names(selection_data))) {
      stopf("Selection z_GLMM.corrected column missing (chr %d, window %s).", num, window_ID)
    }
    sel_z <- suppressWarnings(as.numeric(selection_data[[col_Z_sel]]))
    sel_central <- selection_data[
      POS >= start_central & POS < end_central &
      is.finite(sel_z) & (sel_z^2) > chi_sq_threshold
    ]
    if (nrow(sel_central) == 0L) {
      cat(sprintf("\r[CHR %d] Window %d/%d: no SELECTION hits in central region by z_GLMM.corrected (skip)", num, i, nrow(windows_dt)))
      next
    }

    # GWAS: prefer Z; if absent, fall back to P_value
    central_data <- data.table()
    if (!is.null(col_Z) && col_Z %in% names(gwas_window)) {
      zvec <- suppressWarnings(as.numeric(gwas_window[[col_Z]]))
      central_mask <- (as.integer(gwas_window[[col_POS]]) >= start_central) &
                      (as.integer(gwas_window[[col_POS]]) <  end_central) &
                      (is.finite(zvec) & (zvec^2) > chi_sq_threshold)
      central_data <- gwas_window[central_mask, , drop = FALSE]
    } else if (!is.null(col_PVAL) && col_PVAL %in% names(gwas_window)) {
      pvec <- suppressWarnings(as.numeric(gwas_window[[col_PVAL]]))
      central_mask <- (as.integer(gwas_window[[col_POS]]) >= start_central) &
                      (as.integer(gwas_window[[col_POS]]) <  end_central) &
                      (is.finite(pvec) & pvec <= p_value_threshold & pvec > 0)
      central_data <- gwas_window[central_mask, , drop = FALSE]
    } else {
      stopf("GWAS has neither Z nor P-value columns (chr %d, window %s).", num, window_ID)
    }
    if (nrow(central_data) == 0L) {
      cat(sprintf("\r[CHR %d] Window %d/%d: no GWAS hits in central region by %s (skip)",
                  num, i, nrow(windows_dt),
                  if (!is.null(col_Z) && col_Z %in% names(gwas_window)) "Z" else "P"))
      next
    }

    # --------------------------
    # Prepare Data_Trait2 using EXACT GWAS logic (8 cases + prop-aware QC)
    # --------------------------

    # Deduplicate IDs (keep only IDs that appear once in the window)
    if (is.null(col_ID)) stopf("No SNP ID column found (chr %d, window %s).", num, window_ID)
    idx <- gwas_window[, .I[.N == 1L], by = .(id = get(col_ID))]$V1
    gwas_window <- gwas_window[idx]

    # ---- Prop_cases classification (always classify when column exists) ----
    prop_tag <- "NA"
    sdY_med  <- NA_real_
    s_med    <- NA_real_
    if (!is.null(col_Prop) && (col_Prop %in% names(gwas_window))) {
      pc <- trimws(as.character(gwas_window[[col_Prop]]))

      is_na     <- (pc == "") | toupper(pc) == "NA"
      is_na_sdY <- grepl("^NA_?SDY$", pc, ignore.case = TRUE)
      is_sd     <- grepl("^SD_", pc, ignore.case = TRUE)

      sd_vals <- suppressWarnings(as.numeric(sub("(?i)^SD_\\s*", "", pc[is_sd])))
      leftovers <- pc[!(is_sd | is_na_sdY | is_na)]
      num_vals  <- suppressWarnings(as.numeric(leftovers))

      if (any(!is.na(sd_vals))) {
        prop_tag <- "SD"
        sdY_med  <- median(sd_vals, na.rm = TRUE); if (is.nan(sdY_med)) sdY_med <- NA_real_
      } else if (any(is_na_sdY, na.rm = TRUE)) {
        prop_tag <- "NA_sdY"
      } else if (any(!is.na(num_vals))) {
        prop_tag <- "NUM"
        s_med    <- median(num_vals, na.rm = TRUE); if (is.nan(s_med)) s_med <- NA_real_
      } else {
        prop_tag <- "NA"
      }
    }

    # Determine BETA/SE presence (any finite BETA)
    has_BETA <- !is.null(col_BETA) && (col_BETA %in% names(gwas_window)) &&
                any(is.finite(suppressWarnings(as.numeric(gwas_window[[col_BETA]]))))
    has_SE   <- !is.null(col_SE)   && (col_SE   %in% names(gwas_window))

    # Decide case & required fields (same as fine-mapping; case 6 needs MAF & N)
    case_id <- NA_integer_
    needs <- list(beta=FALSE, se=FALSE, p=FALSE, maf=FALSE, n=FALSE, sdY=FALSE, s=FALSE)
    if (has_BETA && has_SE) {
      if      (prop_tag == "NA")     { case_id <- 1; needs$beta <- TRUE; needs$se <- TRUE; needs$maf <- TRUE; needs$n <- TRUE }
      else if (prop_tag == "NA_sdY") { case_id <- 2; needs$beta <- TRUE; needs$se <- TRUE; needs$sdY <- TRUE }
      else if (prop_tag == "SD")     { case_id <- 3; needs$beta <- TRUE; needs$se <- TRUE; needs$sdY <- TRUE }
      else                           { case_id <- 4; needs$beta <- TRUE; needs$se <- TRUE }
    } else {
      if (is.null(col_PVAL) || !(col_PVAL %in% names(gwas_window))) {
        stopf("No BETA/SE and no P-value column available (chr %d, window %s).", num, window_ID)
      }
      if      (prop_tag == "NA")     { case_id <- 5; needs$p <- TRUE; needs$maf <- TRUE; needs$n <- TRUE }
      else if (prop_tag == "NA_sdY") { case_id <- 6; needs$p <- TRUE; needs$maf <- TRUE; needs$n <- TRUE; needs$sdY <- TRUE }
      else if (prop_tag == "SD")     { case_id <- 7; needs$p <- TRUE; needs$maf <- TRUE; needs$n <- TRUE; needs$sdY <- TRUE }
      else                           { case_id <- 8; needs$p <- TRUE; needs$maf <- TRUE; needs$n <- TRUE; needs$s <- TRUE }
    }

    # Window-level median N (fallback to chromosome median)
    medN_win <- medN_chr
    if (!is.null(col_N) && col_N %in% names(gwas_window)) {
      Nvals_win <- suppressWarnings(as.numeric(gwas_window[[col_N]]))
      medN_win <- median(Nvals_win, na.rm = TRUE)
      if (is.nan(medN_win)) medN_win <- NA_real_
    }

    # ---- Prop-aware QC: filter ONLY columns needed for chosen case ----
    DT <- copy(gwas_window)
    n0 <- nrow(DT)
    keep <- rep(TRUE, n0)

    # P-values needed?
    if (needs$p) {
      if (is.null(col_PVAL) || !(col_PVAL %in% names(DT))) {
        stopf("Case %d requires P-values but column missing (chr %d, window %s).", case_id, num, window_ID)
      }
      if (!is.double(DT[[col_PVAL]])) DT[, (col_PVAL) := as.numeric(get(col_PVAL))]
      tiny <- .Machine$double.xmin * .Machine$double.eps
      DT[is.finite(get(col_PVAL)) & get(col_PVAL) == 0, (col_PVAL) := tiny]
      keep <- keep & (is.finite(DT[[col_PVAL]]) & DT[[col_PVAL]] > 0 & DT[[col_PVAL]] <= 1)
    }

    # Beta/SE needed?
    if (needs$beta || needs$se) {
      if (is.null(col_BETA) || is.null(col_SE) || !(col_BETA %in% names(DT)) || !(col_SE %in% names(DT))) {
        stopf("Case %d requires BETA/SE but column missing (chr %d, window %s).", case_id, num, window_ID)
      }
      DT[, (col_BETA) := as.numeric(get(col_BETA))]
      DT[, (col_SE)   := as.numeric(get(col_SE))]
      keep <- keep & (is.finite(DT[[col_BETA]]) & is.finite(DT[[col_SE]]) & DT[[col_SE]] > 0)
    }

    # MAF needed?
    if (needs$maf) {
      if (is.null(col_MAF) || !(col_MAF %in% names(DT))) {
        stopf("Case %d requires MAF but column missing (chr %d, window %s).", case_id, num, window_ID)
      }
      DT[, (col_MAF) := as.numeric(get(col_MAF))]
      keep <- keep & (is.finite(DT[[col_MAF]]) & DT[[col_MAF]] > 0 & DT[[col_MAF]] < 1)
    }

    DT <- DT[keep]
    if (nrow(DT) == 0L) {
      stopf("After required QC for case %d, no usable variants remain (chr %d, window %s).", case_id, num, window_ID)
    }

    # Always-required ID/POS present?
    if (is.null(col_ID) || !(col_ID %in% names(DT))) stopf("SNP ID column missing after QC (chr %d, window %s).", num, window_ID)
    if (is.null(col_POS) || !(col_POS %in% names(DT))) stopf("POS column missing after QC (chr %d, window %s).", num, window_ID)

    # Scalars needed?
    if (needs$sdY) {
      if (case_id %in% c(3,7) && is.na(sdY_med)) stopf("Case %d requires sdY (from SD_{num}) but median is NA (chr %d, window %s).", case_id, num, window_ID)
      # case 2/6 => sdY = 1
    }
    if (needs$s && is.na(s_med)) stopf("Case %d requires s (median numeric Prop_cases) but it is NA (chr %d, window %s).", case_id, num, window_ID)
    if (needs$n && is.na(medN_win)) stopf("Case %d requires N but median N is NA for this window (chr %d, window %s).", case_id, num, window_ID)

    # Build common vectors
    snp_vec <- as.character(if (!is.null(col_ID)) DT[[col_ID]] else paste0(DT[[col_CHROM]], "_", DT[[col_POS]]))
    pos_vec <- as.integer(DT[[col_POS]])

    # Compose Data_Trait2 per case (1..8)
    if (case_id %in% c(1,2,3,4)) {
      beta_vec    <- as.numeric(DT[[col_BETA]])
      varbeta_vec <- as.numeric(DT[[col_SE]])^2
      if (case_id == 1) {
        maf_vec_win <- as.numeric(DT[[col_MAF]])
        Data_Trait2 <- list(beta=beta_vec, varbeta=varbeta_vec, snp=snp_vec, position=pos_vec,
                            type="quant", MAF=maf_vec_win, N=as.numeric(medN_win))
      } else if (case_id == 2) {
        Data_Trait2 <- list(beta=beta_vec, varbeta=varbeta_vec, snp=snp_vec, position=pos_vec,
                            type="quant", sdY=1)
      } else if (case_id == 3) {
        Data_Trait2 <- list(beta=beta_vec, varbeta=varbeta_vec, snp=snp_vec, position=pos_vec,
                            type="quant", sdY=as.numeric(sdY_med))
      } else { # 4
        Data_Trait2 <- list(beta=beta_vec, varbeta=varbeta_vec, snp=snp_vec, position=pos_vec,
                            type="cc")
      }
    } else {
      pval_vec    <- as.numeric(DT[[col_PVAL]])
      maf_vec_win <- as.numeric(DT[[col_MAF]])
      if (case_id == 5) {
        Data_Trait2 <- list(pvalues=pval_vec, snp=snp_vec, position=pos_vec,
                            type="quant", MAF=maf_vec_win, N=as.numeric(medN_win))
      } else if (case_id == 6) { # NA_sdY, needs MAF+N, sdY=1
        Data_Trait2 <- list(pvalues=pval_vec, snp=snp_vec, position=pos_vec,
                            type="quant", MAF=maf_vec_win, N=as.numeric(medN_win), sdY=1)
      } else if (case_id == 7) {
        Data_Trait2 <- list(pvalues=pval_vec, snp=snp_vec, position=pos_vec,
                            type="quant", MAF=maf_vec_win, N=as.numeric(medN_win), sdY=as.numeric(sdY_med))
      } else { # 8
        Data_Trait2 <- list(pvalues=pval_vec, snp=snp_vec, position=pos_vec,
                            type="cc", MAF=maf_vec_win, N=as.numeric(medN_win), s=as.numeric(s_med))
      }
    }

    # --------------------------
    # Assemble Data_Selection (unchanged)
    # --------------------------
    Data_Selection <- list(
      beta     = as.numeric(selection_data$beta),
      varbeta  = as.numeric((selection_data$se)^2),
      snp      = as.character(selection_data$ID),   # CHR_POS_REF_ALT
      position = as.integer(selection_data$POS),
      type     = "quant",
      sdY      = 1
    )

    # ---- Run coloc (per window) ----
    my.res_2 <- coloc.abf(dataset1 = Data_Selection, dataset2 = Data_Trait2, p1 = 1e-3, p12 = 5e-5)

    # ---- Per-window posterior summaries (no filtering) ----
    summary_dt <- as.data.table(as.list(my.res_2$summary))
    if ("nsnps" %in% names(summary_dt)) summary_dt[, nsnps := NULL]
    summary_dt[, window := window_ID]
    # To preserve original shape expecting gene_id, set gene_id := trait2_name
    summary_dt[, gene_id := trait2_name]
    setnames(summary_dt,
      old = c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","window"),
      new = c("PP_H0","PP_H1","PP_H2","PP_H3_different_causal","PP_colocalization","window"),
      skip_absent = TRUE
    )
    setcolorder(summary_dt, c("window","gene_id","PP_H0","PP_H1","PP_H2","PP_H3_different_causal","PP_colocalization"))
    pp_by_window_list[[length(pp_by_window_list) + 1]] <- summary_dt

    # ---- SNP-level table ordered by PP.H4 (shared causal) ----
    res_snp <- as.data.table(my.res_2$results[order(-my.res_2$results$SNP.PP.H4), ])
    # map POS via selection or trait2 (here trait2 comes from gwas_window)
    map_pos <- unique(rbind(
      selection_data[, .(ID_hg19 = ID, POS)],
      # use gwas_window ID if exists, else compose ID from CHROM_POS
      data.table(ID_hg19 = if (!is.null(col_ID)) as.character(gwas_window[[col_ID]]) else as.character(paste0(gwas_window[[col_CHROM]],"_",gwas_window[[col_POS]])),
                 POS = as.integer(gwas_window[[col_POS]]))
    ))
    res_snp <- merge(res_snp, map_pos, by.x = "snp", by.y = "ID_hg19", all.x = TRUE)

    # window info & rename
    res_snp[, window := window_ID]
    res_snp[, gene_id := trait2_name]
    res_snp[, start_central := start_central]
    res_snp[, end_central := end_central]
    setnames(res_snp, c("snp","SNP.PP.H4"), c("ID_hg19","PIP_shared_causal"), skip_absent = TRUE)
    setcolorder(res_snp, c("window","gene_id","ID_hg19","POS","PIP_shared_causal","start_central","end_central"))

    # store (no postproc)
    snp_by_window_list[[length(snp_by_window_list) + 1]] <- res_snp

    cat(sprintf("\r[CHR %d] Processed window %d/%d", num, i, nrow(windows_dt)))
  } # end windows loop

  cat("\n")

  # ---- Combine (if any) ----
  if (length(snp_by_window_list) == 0L) { cat(sprintf("[CHR %d] No coloc results for any window.\n", num)); next }

  snp_by_window  <- rbindlist(snp_by_window_list, use.names = TRUE, fill = TRUE)
  pp_by_window   <- rbindlist(pp_by_window_list,  use.names = TRUE, fill = TRUE)

  # ------------------ POST-PROCESSING (YOUR RULE) ------------------
  stopifnot(all(c("window","gene_id","ID_hg19","POS","PIP_shared_causal","start_central","end_central") %in% names(snp_by_window)))

  final_snps <- snp_by_window[
    , {
        in_central <- (POS >= start_central & POS < end_central)
        if (any(in_central, na.rm = TRUE)) {
          .SD[in_central][1]
        } else {
          .SD[which.min(PIP_shared_causal)]
        }
      },
    by = .(gene_id, ID_hg19)
  ]

  # Keep only windows whose own lead SNP (variant-level) survives and is still assigned there
  lead_snp_per_window <- function(dt) {
    dt[, .SD[order(-PIP_shared_causal, POS, ID_hg19)][1], by = window][ , .(window, lead_ID = ID_hg19, lead_PIP = PIP_shared_causal)]
  }
  lead_by_window <- lead_snp_per_window(snp_by_window)
  survivors      <- unique(final_snps[, .(ID_hg19, window)])
  kept_windows   <- merge(
    lead_by_window, survivors,
    by.x = c("lead_ID","window"), by.y = c("ID_hg19","window"),
    all.x = TRUE
  )[!is.na(lead_ID), unique(window)]

  pp_by_window_kept <- pp_by_window[window %in% kept_windows]
  final_snps_kept   <- final_snps[window %in% kept_windows]

  setcolorder(final_snps,      names(snp_by_window))
  setcolorder(final_snps_kept, names(snp_by_window))
  final_snps      <- final_snps[order(as.numeric(gsub(".*_", "", window)), -PIP_shared_causal, POS, gene_id, ID_hg19)]
  final_snps_kept <- final_snps_kept[order(as.numeric(gsub(".*_", "", window)), -PIP_shared_causal, POS, gene_id, ID_hg19)]

  # ------------------ SAVE FOUR OUTPUTS ------------------
  outfile_pp_all   <- sprintf("%s/colocalization_results_by_window_chr%d.tsv", outdir, num)
  outfile_pp_final <- sprintf("%s/colocalization_results_final_chr%d.tsv", outdir, num)
  outfile_snp_all  <- sprintf("%s/PP_shared_causal_results_by_window_chr%d.tsv", outdir, num)
  outfile_snp_final<- sprintf("%s/PP_shared_causal_results_final_chr%d.tsv", outdir, num)

  fwrite(pp_by_window, outfile_pp_all, sep = "\t")
  fwrite(pp_by_window_kept, outfile_pp_final, sep = "\t")
  fwrite(snp_by_window, outfile_snp_all, sep = "\t")
  fwrite(final_snps_kept, outfile_snp_final, sep = "\t")

  cat(sprintf("[CHR %d] Completed. Wrote outputs to %s\n", num, outdir))
}

cat("\nAll chromosomes processed successfully!\n")
