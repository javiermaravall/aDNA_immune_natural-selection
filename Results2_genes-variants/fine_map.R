#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(coloc)

# Retrieve target chr_num
args <- commandArgs(trailingOnly = TRUE)  
chr_num  <- as.numeric(args[1])               # convert the first argument to numeric

# Define the threshold values
p_value_threshold <- 5e-8
chi_sq_threshold <- qchisq(1 - p_value_threshold, df = 1)

outdir <- sprintf(
  "fine-mapping/Selection/results"
)

# Create the directory (and all parents if needed)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Paths to global files for this sim
windows_file <- sprintf("fine-mapping/windows/windows.tsv")
selection_file <- sprintf("LMM_withZ_GLMM.corrected.tsv")

# Read windows once (global)
windows_data_all <- fread(windows_file)

# Loop through chromosomes 1-22
for (num in c(chr_num)) {
  cat(sprintf("\nProcessing chromosome %d...\n", num))

  # Filter windows to this chromosome (parse "chr<num>_start_end")
  windows_data <- copy(windows_data_all)
  windows_data[, chr := as.integer(sub("^chr([0-9]+)_.*$", "\\1", window_ID))]
  windows_data <- windows_data[chr == num]
  windows_data[, chr := NULL]

  if (nrow(windows_data) == 0) {
    cat(sprintf("No windows for chromosome %d\n", num))
    next
  }

  # Initialize a list to store results for each window
  results_list <- list()

  # Iterate over each window
  for (i in 1:nrow(windows_data)) {
  #for (i in c(19)) {
    window <- windows_data[i]
    start <- window$start
    end <- window$end
    start_central <- window$start_central_1Mb
    end_central <- window$end_central_1Mb
    window_ID <- window$window_ID

    # Load the selection data for this CHROM & POS range using awk (keep header with NR==1)
    # $1 = CHROM, $2 = POS (tab-delimited)
    cmd <- sprintf("awk -F '\\t' 'NR==1 || ($1==%d && $2>=%d && $2<%d)' %s", num, start, end, shQuote(selection_file))
    selection_data <- fread(cmd = cmd)

    # Skip if no data found for this window
    if (nrow(selection_data) == 0) {
      cat(sprintf("\rProcessed window %d/%d (empty)", i, nrow(windows_data)))
      next
    }

    # Filter for rows meeting the threshold condition within the central 1Mb region
    central_data <- selection_data[POS >= start_central & POS < end_central & (z_GLMM.corrected^2) > chi_sq_threshold]
    if (nrow(central_data) == 0) {
      cat(sprintf("\rProcessed window %d/%d (no significant variants)", i, nrow(windows_data)))
      next
    }

    # Drop rows with missing beta.glmm
    #selection_data <- selection_data[!is.na(beta)]

    # Prepare Data_Selection for fine-mapping
    Data_Selection <- list(
      beta = as.numeric(selection_data$beta),
      varbeta = as.numeric((selection_data$se)^2),
      snp = as.character(selection_data$ID),
      position = as.integer(selection_data$POS),
      type = "quant",
      sdY = 1
    )

    # Run fine-mapping
    my.res <- finemap.abf(dataset = Data_Selection, p1=1e-3)
    
    #print(head(my.res))
    #print(tail(my.res, 1))
    
    # find null row(s)
    ix_null <- which(my.res$snp == "null")
    # replace with "null_{window_ID}"
    my.res$snp[ix_null] <- paste0("null_", window_ID)
    
    #print(tail(my.res, 1))

    my.res <- my.res[order(-my.res$SNP.PP), ]
    
    #print(head(my.res))

    # Skip if no results after filtering
    if (nrow(my.res) == 0) {
      cat(sprintf("\rProcessed window %d/%d (no significant results)", i, nrow(windows_data)))
      next
    }

    # Convert my.res to a data.table
    my.res <- as.data.table(my.res)

    # Map back POS via join
    mapping_data <- selection_data[, .(ID, POS)]
    my.res <- merge(my.res, mapping_data, by.x = "snp", by.y = "ID", all.x = TRUE)

    # Add window info
    my.res[, window := window_ID]
    my.res[, start_central := start_central]
    my.res[, end_central := end_central]

    # Rename/reorder
    setnames(my.res, c("snp", "SNP.PP"), c("ID", "PIP_coloc"))
    setcolorder(my.res, c("window", "ID", "PIP_coloc", "start_central", "end_central", "POS"))

    # Sort by PIP_coloc within window
    setorder(my.res, -PIP_coloc)

    # Store
    results_list[[length(results_list) + 1]] <- my.res

    # Progress
    cat(sprintf("\rProcessed window %d/%d", i, nrow(windows_data)))
  }
  cat("\n")

  # Skip chromosome if no results were found
  if (length(results_list) == 0) {
    cat(sprintf("No significant results found for chromosome %d\n", num))
    next
  }

  # Combine all window results
  combined_results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

  # Post-processing to resolve duplicate IDs across windows
  final_results <- combined_results[, {
    in_any_central <- any(POS >= start_central & POS < end_central, na.rm = TRUE)
    if (in_any_central) {
      central_rows <- .SD[POS >= start_central & POS < end_central]
      central_rows[which.max(PIP_coloc)]
    } else {
      .SD[which.min(PIP_coloc)]
    }
  }, by = ID]

  # Align columns and ordering
  setcolorder(final_results, names(combined_results))
  final_results <- final_results[order(as.numeric(gsub(".*_", "", window)), -PIP_coloc), ]

  # Save outputs
  output_file_by_window <- sprintf("%s/fine-mapping_results_by_window_chr%d.tsv", outdir, num)
  fwrite(combined_results, output_file_by_window, sep = "\t")

  output_file_final <- sprintf("%s/fine-mapping_results_final_chr%d.tsv", outdir, num)
  fwrite(final_results, output_file_final, sep = "\t")

  cat(sprintf("Completed chromosome %d\n", num))
}

cat("\nAll chromosomes processed successfully!\n")
