## Repository overview

This repository contains custom scripts implementing analyses from the manuscript:

**“Ancient DNA reveals that natural selection has upregulated the immune system over the last 10,000 years”**

### Contents

#### 1. `Results1_GWAS`

Scripts for GWAS-level analyses.

- **`functional-enrichment_fine-mapped_variants.py`**  
  Computes functional enrichment of fine-mapped selection variants using baselineLD v2.2 binary annotations.

- **`meta-analysis_S-LDSC_GWAS-traits.py`**  
  Performs random-effects meta-analysis of S-LDSC results across a set of GWAS traits.

#### 2. `Results2_genes-variants`

Scripts for fine-mapping and colocalization analyses.

- **`fine_map.R`**  
  Performs single-causal-variant fine-mapping with `coloc` under a custom windowing protocol.

- **`colocalization.R`**  
  Performs colocalization with `coloc` under a custom windowing protocol.
