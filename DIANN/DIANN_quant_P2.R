#!/usr/bin/env Rscript
# Quick DIA-NN MaxLFQ Extraction
# Simple script to extract MaxLFQ values from filtered DIA-NN output

# 安装缺失的包
install.packages("tibble")
library(tibble)

# 然后重新运行
maxlfq_data <- extract_maxlfq("filtered_data.parquet")

if (!require(pacman)) install.packages("pacman")
pacman::p_load(arrow, dplyr, tidyr, readr, stringr)

extract_maxlfq <- function(input_file) {
  cat("Quick DIA-NN MaxLFQ Extraction\n")
  cat(strrep("-", 50), "\n")
  
  # Load data
  cat("Loading:", input_file, "\n")
  df <- read_parquet(input_file)
  
  cat("Loaded", formatC(nrow(df), format = "d", big.mark = ","), "rows\n")
  cat("Proteins:", formatC(length(unique(df$Protein.Group)), format = "d", big.mark = ","), "\n")
  cat("Samples:", length(unique(df$Run)), "\n")
  
  # Create sample mapping
  sample_mapping <- df %>%
    distinct(Run) %>%
    mutate(Sample_Name = str_remove_all(basename(Run), "\\.(raw|mzML|d)$") %>%
             str_remove_all("_(DIA|dia)$")) %>%
    deframe()
  
  cat("\nExtracting PG.MaxLFQ values...\n")
  
  if (!"PG.MaxLFQ" %in% names(df)) {
    stop("ERROR: PG.MaxLFQ column not found!")
  }
  
  # Pivot to matrix
  result <- df %>%
    select(Protein.Group, Run, PG.MaxLFQ) %>%
    filter(!is.na(PG.MaxLFQ), PG.MaxLFQ > 0) %>%
    pivot_wider(names_from = Run, values_from = PG.MaxLFQ, values_fn = first) %>%
    column_to_rownames("Protein.Group")
  
  # Rename columns
  names(result) <- sample_mapping[names(result)]
  
  # Add annotations
  annotations <- df %>%
    select(Protein.Group, Protein.Names, Genes) %>%
    distinct()
  
  final <- annotations %>%
    left_join(result %>% rownames_to_column("Protein.Group"), by = "Protein.Group") %>%
    column_to_rownames("Protein.Group")
  
  # Calculate statistics
  numeric_cols <- final %>% select_if(is.numeric) %>% names()
  final$N_Samples <- rowSums(!is.na(final[numeric_cols]))
  final$Mean_Intensity <- rowMeans(final[numeric_cols], na.rm = TRUE)
  final$CV_percent <- apply(final[numeric_cols], 1, function(x) {
    if (sum(!is.na(x)) < 2) return(NA)
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
  })
  
  # Save output
  output_file <- str_replace(input_file, "\\.parquet$", "_MaxLFQ.csv")
  write_csv(final, output_file)
  
  cat("\nSaved to:", output_file, "\n")
  cat("Done!\n")
  return(final)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat("Usage: Rscript DIANN_quant_P2.R <filtered_parquet_file>\n")
    quit(status = 1)
  }
  extract_maxlfq(args[1])
}

if (!interactive()) {
  main()
}