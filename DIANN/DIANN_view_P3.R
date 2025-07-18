#!/usr/bin/env Rscript
# DIA-NN Proteomics Analysis Pipeline
# Statistical analysis, visualization, and GO enrichment preparation

if (!require(pacman)) install.packages("pacman")
pacman::p_load(dplyr, readr, ggplot2, pheatmap, corrplot, gridExtra, 
               stringr, tidyr, R6)

DIANNAnalyzer <- R6Class("DIANNAnalyzer",
                         public = list(
                           quant_file = NULL,
                           metadata_file = NULL,
                           output_dir = NULL,
                           proteins = NULL,
                           metadata = NULL,
                           sample_cols = NULL,
                           sample_to_group = NULL,
                           groups = NULL,
                           comparisons = NULL,
                           
                           initialize = function(quant_file, metadata_file) {
                             self$quant_file <- quant_file
                             self$metadata_file <- metadata_file
                             self$output_dir <- file.path(dirname(quant_file), "analysis_results")
                             dir.create(self$output_dir, showWarnings = FALSE, recursive = TRUE)
                             
                             cat(strrep("=", 80), "\n")
                             cat("DIA-NN PROTEOMICS ANALYSIS PIPELINE\n")
                             cat(strrep("=", 80), "\n")
                           },
                           
                           load_data = function() {
                             cat("\n1. Loading data...\n")
                             self$proteins <- read_csv(self$quant_file, show_col_types = FALSE)
                             if ("...1" %in% names(self$proteins)) {
                               self$proteins <- self$proteins %>% column_to_rownames("...1")
                             }
                             
                             self$metadata <- read_csv(self$metadata_file, show_col_types = FALSE)
                             
                             self$sample_cols <- names(self$proteins)[!names(self$proteins) %in% 
                                                                        c("Protein.Names", "Genes", "N_Samples", 
                                                                          "Mean_Intensity", "CV_percent")]
                             
                             self$match_samples()
                           },
                           
                           match_samples = function() {
                             self$sample_to_group <- list()
                             for (i in 1:nrow(self$metadata)) {
                               sample_name <- self$metadata$Sample_Name[i]
                               group <- self$metadata$Group[i]
                               
                               for (col in self$sample_cols) {
                                 if (col == sample_name || grepl(sample_name, col, fixed = TRUE)) {
                                   self$sample_to_group[[col]] <- group
                                   break
                                 }
                               }
                             }
                             
                             self$groups <- sort(unique(unlist(self$sample_to_group)))
                             cat("Groups found:", paste(self$groups, collapse = ", "), "\n")
                           },
                           
                           calculate_statistics = function() {
                             cat("\n3. Calculating statistics...\n")
                             self$comparisons <- list()
                             
                             group_pairs <- combn(self$groups, 2, simplify = FALSE)
                             for (pair in group_pairs) {
                               group1 <- pair[1]
                               group2 <- pair[2]
                               
                               samples1 <- names(self$sample_to_group)[sapply(self$sample_to_group, function(x) x == group1)]
                               samples2 <- names(self$sample_to_group)[sapply(self$sample_to_group, function(x) x == group2)]
                               
                               results <- data.frame()
                               for (protein_id in rownames(self$proteins)) {
                                 values1 <- as.numeric(self$proteins[protein_id, samples1])
                                 values2 <- as.numeric(self$proteins[protein_id, samples2])
                                 
                                 values1 <- values1[!is.na(values1)]
                                 values2 <- values2[!is.na(values2)]
                                 
                                 if (length(values1) < 2 || length(values2) < 2) next
                                 
                                 mean1 <- mean(values1)
                                 mean2 <- mean(values2)
                                 log2fc <- log2(pmax(mean1, 1) / pmax(mean2, 1))
                                 
                                 test_result <- t.test(values1, values2, var.equal = FALSE)
                                 
                                 results <- rbind(results, data.frame(
                                   Protein = protein_id,
                                   Gene = ifelse("Genes" %in% names(self$proteins), 
                                                 self$proteins[protein_id, "Genes"], ""),
                                   Log2FC = log2fc,
                                   P.Value = test_result$p.value,
                                   stringsAsFactors = FALSE
                                 ))
                               }
                               
                               if (nrow(results) > 0) {
                                 results$Q.Value <- p.adjust(results$P.Value, method = "fdr")
                                 results$minus_log10_P <- -log10(results$P.Value)
                                 
                                 comparison_name <- paste(group1, "vs", group2, sep = "_")
                                 write_csv(results, file.path(self$output_dir, paste0("comparison_", comparison_name, ".csv")))
                                 
                                 self$comparisons[[comparison_name]] <- list(
                                   name = comparison_name,
                                   data = results,
                                   group1 = group1,
                                   group2 = group2
                                 )
                               }
                             }
                           },
                           
                           create_visualizations = function() {
                             cat("\n4. Creating visualizations...\n")
                             
                             # PCA plot
                             pca_data <- self$proteins[self$sample_cols] %>% t() %>% as.data.frame() %>% replace(is.na(.), 0)
                             pca_result <- prcomp(pca_data, scale. = TRUE)
                             
                             plot_data <- data.frame(
                               PC1 = pca_result$x[, 1],
                               PC2 = pca_result$x[, 2],
                               Sample = rownames(pca_result$x),
                               Group = sapply(rownames(pca_result$x), function(x) self$sample_to_group[[x]] %||% "Unknown")
                             )
                             
                             p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group)) +
                               geom_point(size = 4) + labs(title = "PCA of Samples") + theme_bw()
                             
                             ggsave(file.path(self$output_dir, "pca_plot.pdf"), p, width = 10, height = 8)
                             
                             # Volcano plots
                             for (comp in self$comparisons) {
                               data <- comp$data
                               data$significance <- "Not Significant"
                               data$significance[data$Q.Value < 0.05 & data$Log2FC > 1] <- "Up"
                               data$significance[data$Q.Value < 0.05 & data$Log2FC < -1] <- "Down"
                               
                               p <- ggplot(data, aes(x = Log2FC, y = minus_log10_P, color = significance)) +
                                 geom_point() + 
                                 scale_color_manual(values = c("Not Significant" = "gray", "Up" = "red", "Down" = "blue")) +
                                 labs(title = paste("Volcano Plot:", comp$name)) + theme_bw()
                               
                               ggsave(file.path(self$output_dir, paste0("volcano_", comp$name, ".pdf")), 
                                      p, width = 10, height = 8)
                             }
                           },
                           
                           prepare_go_analysis = function() {
                             cat("\n5. Preparing for GO analysis...\n")
                             for (comp in self$comparisons) {
                               data <- comp$data
                               up_genes <- data %>% filter(Q.Value < 0.05, Log2FC > 1, !is.na(Gene)) %>% pull(Gene)
                               down_genes <- data %>% filter(Q.Value < 0.05, Log2FC < -1, !is.na(Gene)) %>% pull(Gene)
                               
                               if (length(up_genes) > 0) {
                                 writeLines(up_genes, file.path(self$output_dir, paste0("genes_up_", comp$name, ".txt")))
                               }
                               if (length(down_genes) > 0) {
                                 writeLines(down_genes, file.path(self$output_dir, paste0("genes_down_", comp$name, ".txt")))
                               }
                             }
                           },
                           
                           run_analysis = function() {
                             self$load_data()
                             self$calculate_statistics()
                             self$create_visualizations()
                             self$prepare_go_analysis()
                             cat("\nAnalysis complete!\n")
                           }
                         )
)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    cat("Usage: Rscript DIANN_view_P3.R <maxlfq_file> <metadata_file>\n")
    quit(status = 1)
  }
  
  analyzer <- DIANNAnalyzer$new(args[1], args[2])
  analyzer$run_analysis()
}

if (!interactive()) {
  main()
}