##############################################################################
# WGCNA Analysis for Mass Spectrometry Proteomics Data
# 8 groups (G2_CAPH, G2_CAPH2, G2_KIF4A, G2_TOP2A,
#           PM_CAPH, PM_CAPH2, PM_KIF4A, PM_TOP2A) with replicates
# Targets 6 co-expression modules and exports protein lists per module
##############################################################################

# ── 0. Install / load packages ──────────────────────────────────────────────
required_pkgs <- c("WGCNA", "tidyverse", "pheatmap", "dynamicTreeCut",
                   "flashClust", "openxlsx")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(new_pkgs, ask = FALSE, update = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

allowWGCNAThreads()          # enable multi-threading inside WGCNA

# ── 1. Read & clean data ────────────────────────────────────────────────────
# >>> EDIT THIS PATH to point to your CSV file <<<
input_file <- "your_mass_spec_data.csv"

raw <- read.csv(input_file, header = TRUE, check.names = FALSE,
                stringsAsFactors = FALSE)

cat("Raw dimensions:", nrow(raw), "rows x", ncol(raw), "columns\n")
head(raw[, 1:8])

# Identify metadata vs abundance columns
# Metadata columns: UniprotID, GeneSymbol, Description, Peptides, MosaicID
meta_cols <- c("UniprotID", "GeneSymbol", "Description", "Peptides", "MosaicID")
# Adjust if your column names differ slightly:
meta_cols <- intersect(meta_cols, colnames(raw))
abundance_cols <- setdiff(colnames(raw), meta_cols)

cat("\nMetadata columns found:", paste(meta_cols, collapse = ", "), "\n")
cat("Abundance columns (", length(abundance_cols), "):",
    paste(abundance_cols, collapse = ", "), "\n")

# Remove normalization / header rows (e.g. rows labelled "NORMAL" or all 1s)
if ("UniprotID" %in% colnames(raw)) {
  normal_rows <- grepl("NORMAL|Normal|normal", raw$UniprotID, ignore.case = TRUE)
  if (any(normal_rows)) {
    cat("Removing", sum(normal_rows), "normalization row(s)\n")
    raw <- raw[!normal_rows, ]
  }
}

# Keep protein annotation for later
protein_info <- raw[, meta_cols, drop = FALSE]

# Build numeric abundance matrix  (proteins = rows, samples = columns)
expr_mat <- as.data.frame(lapply(raw[, abundance_cols], as.numeric))
rownames(expr_mat) <- if ("GeneSymbol" %in% meta_cols) {
  make.unique(raw$GeneSymbol)
} else if ("UniprotID" %in% meta_cols) {
  make.unique(raw$UniprotID)
} else {
  paste0("Protein_", seq_len(nrow(raw)))
}

# ── 2. Filtering & imputation ───────────────────────────────────────────────
# Remove proteins with too many missing values (> 50 % NA across samples)
na_frac <- rowMeans(is.na(expr_mat))
keep <- na_frac <= 0.5
cat("\nFiltering: keeping", sum(keep), "of", length(keep),
    "proteins (≤ 50 % NA)\n")
expr_mat <- expr_mat[keep, ]
protein_info <- protein_info[keep, ]

# Impute remaining NAs with the row (protein) minimum (common in proteomics)
expr_mat <- t(apply(expr_mat, 1, function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE)
  return(x)
}))

# Log2-transform (skip if data already log-scaled)
if (max(expr_mat, na.rm = TRUE) > 100) {
  cat("Data appears non-log-transformed; applying log2(x+1)\n")
  expr_mat <- log2(expr_mat + 1)
}

# ── 3. Prepare WGCNA input ──────────────────────────────────────────────────
# WGCNA expects samples (rows) x genes/proteins (columns)
datExpr <- as.data.frame(t(expr_mat))

# Remove proteins with zero variance
goodGenes <- goodSamplesGenes(datExpr, verbose = 3)
if (!goodGenes$allOK) {
  datExpr <- datExpr[goodGenes$goodSamples, goodGenes$goodGenes]
  protein_info <- protein_info[goodGenes$goodGenes, ]
  cat("Removed bad samples/genes; new dims:", dim(datExpr), "\n")
}

cat("\nWGCNA input matrix:", nrow(datExpr), "samples x",
    ncol(datExpr), "proteins\n")

# ── 4. Build sample trait table ──────────────────────────────────────────────
# Parse group names from sample (column) names
sample_names <- rownames(datExpr)
cat("\nSample names:\n"); print(sample_names)

# Extract group label by stripping the replicate suffix
# Handles patterns like "G2_CAPH_1", "PM_TOP2A_Rep2", etc.
group_labels <- sub("_?[Rr]?[Ee]?[Pp]?\\d+$", "", sample_names)
group_labels <- sub("_$", "", group_labels)            # trim trailing _
cat("\nGroup assignments:\n"); print(table(group_labels))

# One-hot encode groups as numeric traits for module-trait correlation
trait_df <- data.frame(SampleName = sample_names,
                       Group = group_labels,
                       stringsAsFactors = FALSE)
trait_design <- model.matrix(~ 0 + Group, data = trait_df)
colnames(trait_design) <- sub("^Group", "", colnames(trait_design))
rownames(trait_design) <- sample_names

# Add phase (G2 vs PM) and bait as extra traits
trait_df$Phase <- ifelse(grepl("^G2", group_labels), 1, 0)
trait_df$Bait  <- group_labels                   # kept for labelling

# ── 5. Pick soft-thresholding power ─────────────────────────────────────────
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         verbose = 5, networkType = "signed")

# Plot scale-free topology fit
pdf("WGCNA_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R²)",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()
cat("\nSaved: WGCNA_soft_threshold.pdf\n")

# Choose power: use WGCNA suggestion or sensible default
chosen_power <- sft$powerEstimate
if (is.na(chosen_power) || chosen_power < 1) chosen_power <- 12
cat("Chosen soft-threshold power:", chosen_power, "\n")

# ── 6. Network construction & module detection ──────────────────────────────
# Using blockwiseModules for memory-efficient one-step network + module ID.
# mergeCutHeight controls how aggressively modules are merged.
# deepSplit and minModuleSize are tuned to aim for ~6 modules.
net <- blockwiseModules(
  datExpr,
  power             = chosen_power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  deepSplit         = 2,
  mergeCutHeight    = 0.25,       # increase to merge more → fewer modules
  pamRespectsDendro = FALSE,
  numericLabels     = TRUE,
  saveTOMs          = FALSE,
  verbose           = 3
)

# Check how many modules we got
n_mods <- length(unique(net$colors)) - 1   # exclude module 0 (grey/unassigned)
cat("\nModules detected:", n_mods, "\n")

# If we got more than 6 modules, progressively merge until we reach 6
if (n_mods > 6) {
  cat("Merging modules to target 6...\n")
  ME_diss <- 1 - cor(net$MEs)
  ME_tree <- hclust(as.dist(ME_diss), method = "average")

  # Iterate over merge heights to find one yielding ≤ 6 modules
  for (h in seq(0.25, 0.90, by = 0.05)) {
    merge_result <- mergeCloseModules(datExpr, net$colors,
                                      cutHeight = h, verbose = 0)
    n_merged <- length(unique(merge_result$colors)) - 1
    if (n_merged <= 6) {
      net$colors <- merge_result$colors
      net$MEs    <- merge_result$newMEs
      cat("Merged at height", h, "→", n_merged, "modules\n")
      break
    }
  }
}

# Convert numeric labels to colour labels
moduleColors <- labels2colors(net$colors)
cat("\nFinal module sizes:\n")
print(table(moduleColors))

# ── 7. Dendrogram + module colour plot ───────────────────────────────────────
pdf("WGCNA_dendrogram_modules.pdf", width = 14, height = 7)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colours",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colours")
dev.off()
cat("Saved: WGCNA_dendrogram_modules.pdf\n")

# ── 8. Module–trait correlation heatmap ──────────────────────────────────────
MEs <- orderMEs(net$MEs)
# Remove the grey (unassigned) module eigengene if present
MEs <- MEs[, !grepl("ME0|MEgrey", colnames(MEs)), drop = FALSE]

modTraitCor  <- cor(MEs, trait_design, use = "p")
modTraitPval <- corPvalueStudent(modTraitCor, nrow(datExpr))

# Build text matrix with r (p)
textMatrix <- paste0(signif(modTraitCor, 2), "\n(",
                     signif(modTraitPval, 1), ")")
dim(textMatrix) <- dim(modTraitCor)

pdf("WGCNA_module_trait_heatmap.pdf", width = 12, height = 8)
par(mar = c(8, 10, 3, 3))
labeledHeatmap(Matrix    = modTraitCor,
               xLabels   = colnames(trait_design),
               yLabels   = colnames(MEs),
               ySymbols  = colnames(MEs),
               colorLabels = FALSE,
               colors    = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text  = 0.6,
               zlim      = c(-1, 1),
               main      = "Module–trait relationships")
dev.off()
cat("Saved: WGCNA_module_trait_heatmap.pdf\n")

# ── 9. Module membership & hub proteins ──────────────────────────────────────
geneModuleMembership <- cor(datExpr, MEs, use = "p")
MMPvalue             <- corPvalueStudent(geneModuleMembership, nrow(datExpr))

colnames(geneModuleMembership) <- paste0("MM.", gsub("^ME", "", colnames(MEs)))
colnames(MMPvalue)             <- paste0("p.MM.", gsub("^ME", "", colnames(MEs)))

# ── 10. Export protein lists per module ──────────────────────────────────────
output_dir <- "WGCNA_module_protein_lists"
dir.create(output_dir, showWarnings = FALSE)

unique_modules <- sort(unique(moduleColors))
unique_modules <- unique_modules[unique_modules != "grey"]  # exclude unassigned

module_summary <- list()

for (mod in unique_modules) {
  idx <- which(moduleColors == mod)
  mod_proteins <- colnames(datExpr)[idx]

  # Build a data frame with protein info + module membership
  out_df <- data.frame(
    ProteinName = mod_proteins,
    stringsAsFactors = FALSE
  )

  # Try to append annotation columns
  gene_idx <- match(mod_proteins, rownames(expr_mat))
  if (length(meta_cols) > 0 && !all(is.na(gene_idx))) {
    for (mc in meta_cols) {
      out_df[[mc]] <- protein_info[[mc]][gene_idx]
    }
  }

  # Add module membership score for this module
  mm_col <- paste0("MM.", mod)
  if (mm_col %in% colnames(geneModuleMembership)) {
    out_df$ModuleMembership <- geneModuleMembership[idx, mm_col]
    out_df <- out_df[order(-abs(out_df$ModuleMembership)), ]
  }

  fname <- file.path(output_dir, paste0("Module_", mod, "_proteins.csv"))
  write.csv(out_df, fname, row.names = FALSE)
  cat("Module", mod, ":", nrow(out_df), "proteins →", fname, "\n")

  module_summary[[mod]] <- out_df
}

# Also save all modules together in one Excel workbook
wb <- createWorkbook()
for (mod in names(module_summary)) {
  addWorksheet(wb, mod)
  writeData(wb, mod, module_summary[[mod]])
}
# Add a summary sheet
summary_df <- data.frame(
  Module     = names(module_summary),
  N_Proteins = sapply(module_summary, nrow)
)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
saveWorkbook(wb, "WGCNA_all_modules.xlsx", overwrite = TRUE)
cat("\nSaved: WGCNA_all_modules.xlsx\n")

# ── 11. Module eigengene bar plots across groups ─────────────────────────────
ME_long <- MEs %>%
  mutate(Sample = rownames(MEs),
         Group  = group_labels) %>%
  pivot_longer(cols = starts_with("ME"),
               names_to = "Module", values_to = "Eigengene")

pdf("WGCNA_eigengene_barplots.pdf", width = 14, height = 10)
ggplot(ME_long, aes(x = Group, y = Eigengene, fill = Group)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module eigengene expression across groups",
       y = "Module eigengene value") +
  scale_fill_brewer(palette = "Set2")
dev.off()
cat("Saved: WGCNA_eigengene_barplots.pdf\n")

# ── 12. Summary ──────────────────────────────────────────────────────────────
cat("\n",
    "═══════════════════════════════════════════════════════\n",
    "  WGCNA analysis complete\n",
    "═══════════════════════════════════════════════════════\n",
    "  Soft-threshold power : ", chosen_power, "\n",
    "  Total proteins used  : ", ncol(datExpr), "\n",
    "  Modules detected     : ", length(unique_modules), "\n",
    "  Unassigned (grey)    : ", sum(moduleColors == "grey"), "\n\n",
    "  Output files:\n",
    "   • WGCNA_soft_threshold.pdf        - power selection plots\n",
    "   • WGCNA_dendrogram_modules.pdf    - dendrogram + module colours\n",
    "   • WGCNA_module_trait_heatmap.pdf  - module–trait correlations\n",
    "   • WGCNA_eigengene_barplots.pdf    - eigengenes across groups\n",
    "   • WGCNA_all_modules.xlsx          - all module protein lists\n",
    "   • ", output_dir, "/                - individual module CSVs\n",
    "═══════════════════════════════════════════════════════\n",
    sep = "")
