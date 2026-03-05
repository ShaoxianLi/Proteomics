#!/usr/bin/env Rscript
# ==============================================================================
# Protein Coverage Fraction Distribution Analysis
# ==============================================================================

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cran.r-project.org")
library(ggplot2)

# ---- Read TSV file (modify path as needed) ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  tsv_file <- args[1]
} else {
  tsv_file <- "your_file.tsv"  # <-- Change this to your TSV file path
}

df <- read.delim(tsv_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# ---- Basic statistics ----
cov <- df$coverage_fraction
cat("===== Coverage Fraction Summary =====\n")
cat(sprintf("Total proteins:  %d\n", length(cov)))
cat(sprintf("Mean:            %.4f (%.2f%%)\n", mean(cov), mean(cov) * 100))
cat(sprintf("Median:          %.4f (%.2f%%)\n", median(cov), median(cov) * 100))
cat(sprintf("SD:              %.4f\n", sd(cov)))
cat(sprintf("Min:             %.4f (%.2f%%)\n", min(cov), min(cov) * 100))
cat(sprintf("Max:             %.4f (%.2f%%)\n", max(cov), max(cov) * 100))
cat(sprintf("Q1 (25%%):        %.4f\n", quantile(cov, 0.25)))
cat(sprintf("Q3 (75%%):        %.4f\n", quantile(cov, 0.75)))
cat("=====================================\n")

mean_val   <- mean(cov)
median_val <- median(cov)

# ---- Plot 1: Histogram with density overlay ----
p1 <- ggplot(df, aes(x = coverage_fraction)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "#377EB8", color = "white", alpha = 0.8) +
  geom_density(linewidth = 0.8, color = "#E41A1C") +
  geom_vline(xintercept = mean_val, linetype = "dashed", color = "#FF7F00", linewidth = 0.9) +
  geom_vline(xintercept = median_val, linetype = "dotted", color = "#4DAF4A", linewidth = 0.9) +
  annotate("text", x = mean_val + 0.02, y = Inf, vjust = 2,
           label = sprintf("Mean = %.3f", mean_val), color = "#FF7F00", size = 4, fontface = "bold") +
  annotate("text", x = median_val - 0.02, y = Inf, vjust = 3.5,
           label = sprintf("Median = %.3f", median_val), color = "#4DAF4A", size = 4, fontface = "bold") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x     = "Coverage Fraction",
    y     = "Density",
    title = "Protein Coverage Fraction Distribution"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ---- Plot 2: Cumulative distribution (ECDF) ----
p2 <- ggplot(df, aes(x = coverage_fraction)) +
  stat_ecdf(geom = "step", linewidth = 0.8, color = "#377EB8") +
  geom_vline(xintercept = mean_val, linetype = "dashed", color = "#FF7F00", linewidth = 0.9) +
  geom_vline(xintercept = median_val, linetype = "dotted", color = "#4DAF4A", linewidth = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50", linewidth = 0.5) +
  annotate("text", x = mean_val + 0.02, y = 0.95,
           label = sprintf("Mean = %.3f", mean_val), color = "#FF7F00", size = 4, fontface = "bold") +
  annotate("text", x = median_val - 0.02, y = 0.45,
           label = sprintf("Median = %.3f", median_val), color = "#4DAF4A", size = 4, fontface = "bold") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = "Coverage Fraction",
    y     = "Cumulative Proportion",
    title = "Cumulative Distribution of Coverage Fraction (ECDF)"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ---- Plot 3: Box + Violin plot ----
p3 <- ggplot(df, aes(x = "", y = coverage_fraction)) +
  geom_violin(fill = "#377EB8", alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.15, fill = "#377EB8", alpha = 0.6, outlier.shape = 21,
               outlier.fill = "#E41A1C", outlier.size = 2) +
  geom_hline(yintercept = mean_val, linetype = "dashed", color = "#FF7F00", linewidth = 0.8) +
  geom_hline(yintercept = median_val, linetype = "dotted", color = "#4DAF4A", linewidth = 0.8) +
  annotate("text", x = 1.35, y = mean_val,
           label = sprintf("Mean = %.3f", mean_val), color = "#FF7F00", size = 4, fontface = "bold") +
  annotate("text", x = 1.35, y = median_val - 0.03,
           label = sprintf("Median = %.3f", median_val), color = "#4DAF4A", size = 4, fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = "",
    y     = "Coverage Fraction",
    title = "Coverage Fraction (Violin + Boxplot)"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ---- Save plots ----
ggsave("coverage_histogram.png", p1, width = 8, height = 5, dpi = 300)
ggsave("coverage_ecdf.png",      p2, width = 8, height = 5, dpi = 300)
ggsave("coverage_violin.png",    p3, width = 5, height = 6, dpi = 300)

ggsave("coverage_histogram.pdf", p1, width = 8, height = 5)
ggsave("coverage_ecdf.pdf",      p2, width = 8, height = 5)
ggsave("coverage_violin.pdf",    p3, width = 5, height = 6)

cat("\nPlots saved: coverage_histogram, coverage_ecdf, coverage_violin (.png & .pdf)\n")
