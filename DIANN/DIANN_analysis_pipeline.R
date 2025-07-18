#!/usr/bin/env Rscript
# DIA-NN Master Pipeline
# Orchestrates the complete DIA-NN analysis workflow

# Source other scripts
source_if_exists <- function(file) {
  if (file.exists(file)) {
    source(file)
    return(TRUE)
  }
  return(FALSE)
}

# Try to source other scripts
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
if (script_dir == "") script_dir <- "."

scripts_loaded <- c(
  source_if_exists(file.path(script_dir, "diann_parquet_filtering_P1.R")),
  source_if_exists(file.path(script_dir, "DIANN_quant_P2.R")),
  source_if_exists(file.path(script_dir, "DIANN_view_P3.R"))
)

if (!all(scripts_loaded)) {
  cat("Warning: Some scripts could not be loaded. Make sure all 4 scripts are in the same directory.\n")
}

run_diann_pipeline <- function(input_parquet, 
                               metadata_file,
                               output_dir = NULL,
                               min_peak_width = 0.05,
                               min_pg_maxlfq = NULL) {
  
  cat(strrep("=", 80), "\n")
  cat("DIA-NN MASTER PIPELINE (R VERSION)\n")
  cat(strrep("=", 80), "\n")
  
  # Set output directory
  if (is.null(output_dir)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path(dirname(input_parquet), paste0("DIANN_Analysis_", timestamp))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Step 1: Filter data
  cat("\nSTEP 1: Filtering parquet file\n")
  filtered_parquet <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(input_parquet)), "_filtered.parquet"))
  filter_diann_parquet(input_parquet, filtered_parquet, min_peak_width = min_peak_width, min_pg_maxlfq = min_pg_maxlfq)
  
  # Step 2: Extract MaxLFQ
  cat("\nSTEP 2: Extracting MaxLFQ values\n")
  extract_maxlfq(filtered_parquet)
  
  # Step 3: Analysis
  cat("\nSTEP 3: Statistical analysis\n")
  maxlfq_csv <- str_replace(filtered_parquet, "\\.parquet$", "_MaxLFQ.csv")
  metadata_copy <- file.path(output_dir, "metadata.csv")
  file.copy(metadata_file, metadata_copy)
  
  analyzer <- DIANNAnalyzer$new(maxlfq_csv, metadata_copy)
  analyzer$run_analysis()
  
  cat("\nPipeline completed! Results in:", output_dir, "\n")
  return(output_dir)
}

# 交互式运行时的参数设置
if (interactive()) {
  # 在这里定义你的输入文件路径
  input_parquet <- "/path/to/your/input_file.parquet"
  metadata_file <- "/path/to/your/metadata.csv"
  output_dir <- "/path/to/your/output_directory"  # 可选
  
  # 运行分析
  run_diann_pipeline(input_parquet, metadata_file, output_dir)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    cat("Usage: Rscript DIANN_analysis_pipeline.R <parquet_file> <metadata_file> [options]\n")
    quit(status = 1)
  }
  
  run_diann_pipeline(args[1], args[2])
}

if (!interactive()) {
  main()
}