#!/usr/bin/env Rscript
# DIA-NN Parquet Filtering Script with Peak Width and Intensity Filters
# Filters DIA-NN output based on quality criteria including peak width and intensity

# Load required libraries
if (!require(pacman)) install.packages("pacman")
pacman::p_load(arrow, dplyr, stringr, readr)

# Define operators
`%R%` <- function(string, times) strrep(string, times)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x

filter_diann_parquet <- function(input_file, 
                                 output_file = NULL,
                                 report_file = NULL,
                                 min_peak_width = 0.05,
                                 min_precursor_quantity = NULL,
                                 min_pg_maxlfq = NULL,
                                 min_pg_quantity = NULL) {
  
  cat("=" %R% 60, "\n")
  cat("DIA-NN FILTERING (DIA-SPECIFIC)\n")
  cat("=" %R% 60, "\n")
  
  # Set default output files
  if (is.null(output_file)) {
    output_file <- str_replace(input_file, "\\.parquet$", "_filtered.parquet")
  }
  if (is.null(report_file)) {
    report_file <- str_replace(input_file, "\\.parquet$", "_filter_report.txt")
  }
  
  cat("Input:", input_file, "\n")
  cat("Output:", output_file, "\n")
  cat("Min Peak Width:", min_peak_width, "minutes\n")
  
  # Load data
  cat("\nLoading parquet file...\n")
  df <- read_parquet(input_file)
  initial_rows <- nrow(df)
  cat("Loaded", formatC(initial_rows, format = "d", big.mark = ","), "rows\n")
  
  # Initialize report
  report <- c(
    "DIA-NN Filtering Report (DIA-Specific)",
    strrep("=", 60),
    paste("Generated:", Sys.time()),
    paste("Input file:", input_file),
    paste("Initial rows:", formatC(initial_rows, format = "d", big.mark = ",")),
    paste("Min peak width threshold:", min_peak_width, "minutes"),
    ""
  )
  
  # Calculate peak width
  if (all(c("RT.Start", "RT.Stop") %in% names(df))) {
    df <- df %>% mutate(Peak.Width = RT.Stop - RT.Start)
    cat("Peak width calculated\n")
    report <- c(report, "Peak width statistics added", "")
  }
  
  # Apply filters
  cat("\nApplying filters:\n")
  df_filtered <- df
  
  # Filter definitions
  filters <- list(
    list(name = "Remove contaminants", col = "Genes", 
         condition = function(x) !str_starts(x %||% "", "##")),
    list(name = "Q.Value <= 0.01", col = "Q.Value", 
         condition = function(x) x <= 0.01),
    list(name = "Lib.Q.Value < 0.01", col = "Lib.Q.Value", 
         condition = function(x) x < 0.01),
    list(name = "PG.Q.Value < 0.05", col = "PG.Q.Value", 
         condition = function(x) x < 0.05),
    list(name = "Quantity.Quality > 0.5", col = "Quantity.Quality", 
         condition = function(x) x > 0.5),
    list(name = "PG.MaxLFQ.Quality > 0.7", col = "PG.MaxLFQ.Quality", 
         condition = function(x) x > 0.7)
  )
  
  # Add conditional filters
  if ("Peak.Width" %in% names(df_filtered)) {
    filters <- append(filters, list(
      list(name = paste("Peak.Width >=", min_peak_width), 
           col = "Peak.Width", 
           condition = function(x) x >= min_peak_width)
    ))
  }
  
  if (!is.null(min_precursor_quantity) && "Precursor.Quantity" %in% names(df_filtered)) {
    filters <- append(filters, list(
      list(name = paste("Precursor.Quantity >", min_precursor_quantity),
           col = "Precursor.Quantity",
           condition = function(x) x > min_precursor_quantity)
    ))
  }
  
  if (!is.null(min_pg_maxlfq) && "PG.MaxLFQ" %in% names(df_filtered)) {
    filters <- append(filters, list(
      list(name = paste("PG.MaxLFQ >", min_pg_maxlfq),
           col = "PG.MaxLFQ",
           condition = function(x) x > min_pg_maxlfq)
    ))
  }
  
  # Apply each filter
  for (filter in filters) {
    if (filter$col %in% names(df_filtered)) {
      rows_before <- nrow(df_filtered)
      if (filter$col == "Genes") {
        mask <- filter$condition(df_filtered[[filter$col]])
      } else {
        mask <- filter$condition(df_filtered[[filter$col]]) & !is.na(df_filtered[[filter$col]])
      }
      df_filtered <- df_filtered[mask, ]
      rows_after <- nrow(df_filtered)
      removed <- rows_before - rows_after
      
      cat("  ", filter$name, ":", removed, "removed\n")
      report <- c(report, paste(filter$name, ":", removed, "removed"))
    }
  }
  
  # Summary
  final_removed <- initial_rows - nrow(df_filtered)
  final_pct <- round(final_removed / initial_rows * 100, 1)
  
  cat("\nTotal removed:", final_removed, paste0("(", final_pct, "%)\n"))
  report <- c(report, "", "SUMMARY",
              paste("Final rows:", nrow(df_filtered)),
              paste("Total removed:", final_removed, paste0("(", final_pct, "%)")))
  
  # Remove Peak.Width before saving
  if ("Peak.Width" %in% names(df_filtered)) {
    df_filtered <- df_filtered %>% select(-Peak.Width)
  }
  
  # Save files
  write_parquet(df_filtered, output_file)
  writeLines(report, report_file)
  
  cat("Filtering complete!\n")
  return(df_filtered)
}

# Main function for command line usage
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript diann_parquet_filtering_P1.R <input_file> [options]\n")
    cat("Options:\n")
    cat("  --output OUTPUT_FILE\n")
    cat("  --min-peak-width WIDTH\n")
    cat("  --min-pg-maxlfq VAL\n")
    quit(status = 1)
  }
  
  # Parse arguments (simplified)
  input_file <- args[1]
  output_file <- NULL
  min_peak_width <- 0.05
  min_pg_maxlfq <- NULL
  
  # Simple argument parsing
  i <- 2
  while (i <= length(args)) {
    if (args[i] == "--output" && i + 1 <= length(args)) {
      output_file <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--min-peak-width" && i + 1 <= length(args)) {
      min_peak_width <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--min-pg-maxlfq" && i + 1 <= length(args)) {
      min_pg_maxlfq <- as.numeric(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  filter_diann_parquet(input_file, output_file, min_peak_width = min_peak_width,
                       min_pg_maxlfq = min_pg_maxlfq)
}

if (!interactive()) {
  main()
}