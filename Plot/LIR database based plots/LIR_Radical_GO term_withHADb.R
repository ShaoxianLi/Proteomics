# ============= 独立运行的完整SL18分析代码 (英文版) =============
# This code can be run independently and generates all plots in English

# ============= Step 1: Install and Load Required Packages =============

required_packages <- c("dplyr", "ggplot2", "tidyr", "RColorBrewer")
bioc_packages <- c("org.Hs.eg.db", "GO.db", "AnnotationDbi")

# Install CRAN packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  cat("Installing CRAN packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

# Install Bioconductor packages
missing_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[,"Package"]]
if(length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_bioc)
}

# Load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)

cat("All packages loaded successfully!\n")

# ============= Step 2: Read Your Real SL18 Data =============

# Your data file path
data_file_path <- "C:/Users/argo_/OneDrive - UMass Chan Medical School/Desktop/Result/plot/SL18_20250625_APEX_LC3B diff autophagy treatment on Eclipse/SL18_Live_pq_153.csv"

cat("Reading SL18 data file...\n")

if (file.exists(data_file_path)) {
  # Read raw data
  raw_data <- read.csv(data_file_path, stringsAsFactors = FALSE)
  cat("Successfully read data, rows:", nrow(raw_data), "\n")
  
  # Remove normalization factor rows
  raw_data <- raw_data[raw_data$MosaicID != "**NORMALIZATION_FACTORS**", ]
  
  # Extract sample column names
  control_cols <- grep("^Ctr_", colnames(raw_data), value = TRUE)
  treatment_cols <- list(
    Torin1_Baf = grep("^Torin1_Baf_", colnames(raw_data), value = TRUE),
    MLSA5_Baf = grep("^MLSA5_Baf_", colnames(raw_data), value = TRUE),
    OA_Baf = grep("^OA_Baf_", colnames(raw_data), value = TRUE),
    Mon_Baf = grep("^Mon_Baf_", colnames(raw_data), value = TRUE),
    Tunica_Baf = grep("^Tunica_Baf_", colnames(raw_data), value = TRUE)
  )
  
  cat("Sample columns identified:\n")
  cat("- Control group:", length(control_cols), "samples\n")
  for (treat in names(treatment_cols)) {
    cat("- ", treat, ":", length(treatment_cols[[treat]]), "samples\n")
  }
  
} else {
  stop("Data file not found, please check path:", data_file_path)
}

# ============= Step 3: Data Processing and Differential Analysis =============

cat("Performing differential protein analysis...\n")

# Calculate fold change for each treatment group vs control
calculate_fold_change <- function(data, control_cols, treatment_cols, treatment_name) {
  # Calculate mean for control and treatment groups
  control_mean <- rowMeans(data[, control_cols], na.rm = TRUE)
  treatment_mean <- rowMeans(data[, treatment_cols], na.rm = TRUE)
  
  # Calculate fold change (log2)
  fold_change <- treatment_mean - control_mean
  
  # Simple t-test
  p_values <- apply(data, 1, function(row) {
    ctrl_vals <- as.numeric(row[control_cols])
    treat_vals <- as.numeric(row[treatment_cols])
    
    # Remove NA values
    ctrl_vals <- ctrl_vals[!is.na(ctrl_vals)]
    treat_vals <- treat_vals[!is.na(treat_vals)]
    
    if (length(ctrl_vals) >= 2 && length(treat_vals) >= 2) {
      t_test <- t.test(treat_vals, ctrl_vals)
      return(t_test$p.value)
    } else {
      return(1)
    }
  })
  
  result <- data.frame(
    GeneSymbol = data$GeneSymbol,
    UniprotID = data$UniprotID,
    Description = data$Description,
    Treatment = treatment_name,
    FoldChange = fold_change,
    LogFoldChange = fold_change,  # Already log2 values
    PValue = p_values,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Calculate fold change for each treatment group
all_comparisons <- data.frame()

for (treat_name in names(treatment_cols)) {
  comparison_result <- calculate_fold_change(
    raw_data, 
    control_cols, 
    treatment_cols[[treat_name]], 
    treat_name
  )
  all_comparisons <- rbind(all_comparisons, comparison_result)
}

# Filter significantly upregulated proteins
significance_threshold <- 0.05
fold_change_threshold <- 0.5  # log2(1.4) ≈ 0.5

all_significant_proteins <- all_comparisons %>%
  dplyr::filter(
    PValue < significance_threshold,
    LogFoldChange > fold_change_threshold
  ) %>%
  dplyr::mutate(
    FoldChange = 2^LogFoldChange  # Convert to linear fold change for visualization
  )

cat("Differential analysis complete!\n")
cat("Significantly upregulated proteins by treatment:\n")
table(all_significant_proteins$Treatment) %>% print()

# ============= Step 4: Get Real GO Annotations =============

cat("Getting real GO annotations from database...\n")

# Get GO annotations for all genes
get_go_annotations <- function(gene_symbols) {
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = gene_symbols,
                       column = "ENTREZID", 
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Remove NA values
  valid_entrez <- entrez_ids[!is.na(entrez_ids)]
  
  if (length(valid_entrez) == 0) {
    return(data.frame())
  }
  
  # Get GO annotations
  go_annotations <- select(org.Hs.eg.db,
                           keys = valid_entrez,
                           columns = c("SYMBOL", "GO", "ONTOLOGY"),
                           keytype = "ENTREZID")
  
  return(go_annotations)
}

# Get GO annotations for all proteins
all_genes <- unique(all_significant_proteins$GeneSymbol)
go_data <- get_go_annotations(all_genes)

cat("GO annotation retrieval complete, annotation entries:", nrow(go_data), "\n")

# Define organelle-related GO terms
er_go_terms <- c(
  "GO:0005783",  # endoplasmic reticulum
  "GO:0005788",  # endoplasmic reticulum lumen
  "GO:0005789",  # endoplasmic reticulum membrane
  "GO:0005793",  # endoplasmic reticulum-Golgi intermediate compartment
  "GO:0030176",  # integral component of endoplasmic reticulum membrane
  "GO:0071556",  # integral component of lumenal side of endoplasmic reticulum membrane
  "GO:0005790",  # smooth endoplasmic reticulum
  "GO:0005791"   # rough endoplasmic reticulum
)

golgi_go_terms <- c(
  "GO:0005794",  # Golgi apparatus
  "GO:0005795",  # Golgi stack
  "GO:0005796",  # Golgi lumen
  "GO:0005797",  # Golgi medial cisterna
  "GO:0005798",  # Golgi-associated vesicle
  "GO:0005799",  # Golgi-associated vesicle membrane
  "GO:0000139",  # Golgi membrane
  "GO:0005801",  # cis-Golgi network
  "GO:0005802"   # trans-Golgi network
)

mito_go_terms <- c(
  "GO:0005739",  # mitochondrion
  "GO:0005740",  # mitochondrial envelope
  "GO:0005741",  # mitochondrial outer membrane
  "GO:0005743",  # mitochondrial inner membrane
  "GO:0005744",  # mitochondrial intermembrane space
  "GO:0005759",  # mitochondrial matrix
  "GO:0005758",  # mitochondrial intermembrane space
  "GO:0005746",  # mitochondrial respirasome
  "GO:0005742"   # mitochondrial outer membrane translocase complex
)

autophagy_go_terms <- c(
  "GO:0006914",  # autophagy
  "GO:0000045",  # autophagosome assembly
  "GO:0034045",  # autophagosome membrane
  "GO:0005776",  # autophagosome
  "GO:0000421",  # autophagosome membrane
  "GO:0016236",  # macroautophagy
  "GO:0061908",  # phagophore
  "GO:0034274",  # autophagosome maturation
  "GO:0097352",  # autophagosome maturation
  "GO:0006995"   # cellular response to nitrogen starvation
)

# Create GO annotation databases
create_go_protein_list <- function(go_terms, go_data) {
  proteins <- go_data %>%
    dplyr::filter(GO %in% go_terms) %>%
    dplyr::distinct(SYMBOL) %>%
    dplyr::pull(SYMBOL)
  
  return(data.frame(SYMBOL = proteins, stringsAsFactors = FALSE))
}

df_er_prots_all <- create_go_protein_list(er_go_terms, go_data)
df_golgi_prots_all <- create_go_protein_list(golgi_go_terms, go_data)
df_mito_prots_all <- create_go_protein_list(mito_go_terms, go_data)
df_autophagy_prots_all <- create_go_protein_list(autophagy_go_terms, go_data)

cat("Real GO annotation databases created:\n")
cat("- ER proteins:", nrow(df_er_prots_all), "entries\n")
cat("- Golgi proteins:", nrow(df_golgi_prots_all), "entries\n")
cat("- Mitochondrial proteins:", nrow(df_mito_prots_all), "entries\n")
cat("- Autophagy proteins:", nrow(df_autophagy_prots_all), "entries\n")

# ============= Step 5: Create LIR Protein Database =============

# Try to read your LIR database file
possible_lir_paths <- c(
  "C:/Users/argo_/OneDrive - UMass Chan Medical School/Desktop/Result/database/LIR_Human with score.csv"
)

lir_file_found <- FALSE
for (path in possible_lir_paths) {
  if (file.exists(path)) {
    lir_file_path <- path
    lir_file_found <- TRUE
    break
  }
}

if (lir_file_found) {
  cat("Reading real LIR database...\n")
  lir_raw_data <- read.csv(lir_file_path, stringsAsFactors = FALSE)
  
  # Process LIR data
  df_lir_proteins <- lir_raw_data %>%
    dplyr::select(UniprotID, GeneSymbol, anchor_no, protein_desc, MF_class, BP_class, CC_class) %>%
    dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
    dplyr::mutate(
      # Create confidence scores based on anchor_no
      Confidence = dplyr::case_when(
        anchor_no >= 15 ~ "High",
        anchor_no >= 10 ~ "Medium", 
        anchor_no >= 5 ~ "Low",
        TRUE ~ "Very_Low"
      ),
      # Create LIR type classification based on GO annotations
      LIR_Type = dplyr::case_when(
        grepl("autophagy|autophag", protein_desc, ignore.case = TRUE) ~ "Autophagy_Core",
        grepl("ubiquitin|ligase", protein_desc, ignore.case = TRUE) ~ "Ubiquitin_System",
        grepl("ribosom|translation", protein_desc, ignore.case = TRUE) ~ "Translation",
        grepl("transcription|DNA", protein_desc, ignore.case = TRUE) ~ "Transcription",
        grepl("kinase|phosphatase", protein_desc, ignore.case = TRUE) ~ "Signaling",
        grepl("membrane|transport", protein_desc, ignore.case = TRUE) ~ "Membrane_Transport",
        grepl("mitochondr", protein_desc, ignore.case = TRUE) ~ "Mitochondrial",
        grepl("cytoskeleton|actin|tubulin", protein_desc, ignore.case = TRUE) ~ "Cytoskeleton",
        grepl("metaboli|enzyme", protein_desc, ignore.case = TRUE) ~ "Metabolism",
        !is.na(protein_desc) & protein_desc != "" ~ "Other_Functional",
        TRUE ~ "Unknown"
      ),
      SYMBOL = GeneSymbol,
      UNIPROT = UniprotID
    ) %>%
    dplyr::filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
    dplyr::arrange(desc(anchor_no))
  
  cat("Real LIR database loaded, protein count:", nrow(df_lir_proteins), "\n")
  
} else {
  cat("LIR file not found, using literature-based LIR database...\n")
  
  # Create literature-based LIR database
  df_lir_proteins <- data.frame(
    SYMBOL = c("SQSTM1", "NBR1", "OPTN", "TAX1BP1", "CALCOCO2", "TOLLIP", "CCPG1",
               "UBR4", "USP5", "BAG3", "DENND4C", "OTUD4", "SPECC1L", "RUFY1", 
               "DNM2", "RABEP2", "PTGES3", "KPNA6", "RABGAP1", "EZR"),
    anchor_no = c(24, 14, 12, 10, 16, 11, 13, 16, 16, 11, 15, 10, 12, 13, 
                  15, 10, 18, 11, 15, 7),
    Confidence = c("High", "Medium", "Medium", "Medium", "High", "Medium", "Medium",
                   "High", "High", "Medium", "High", "Medium", "Medium", "Medium",
                   "High", "Medium", "High", "Medium", "High", "Low"),
    LIR_Type = c(rep("Autophagy_Core", 7), rep("Other_Functional", 13)),
    stringsAsFactors = FALSE
  )
}

# ============= Step 6: Enhanced Protein Annotation with Data Cleaning =============

cat("Annotating proteins with GO and LIR information...\n")

# 首先清理significant proteins数据
all_significant_proteins_clean <- all_significant_proteins %>%
  dplyr::filter(
    !is.na(GeneSymbol),                    # 移除GeneSymbol为NA的行
    GeneSymbol != "",                      # 移除空字符串
    GeneSymbol != "NA",                    # 移除字符串"NA"
    !is.na(UniprotID),                     # 移除UniprotID为NA的行
    !is.na(LogFoldChange),                 # 移除LogFoldChange为NA的行
    !is.infinite(LogFoldChange),           # 移除无限值
    !is.na(PValue)                         # 移除PValue为NA的行
  ) %>%
  dplyr::mutate(
    GeneSymbol = as.character(GeneSymbol),          # 确保是字符类型
    GeneSymbol = trimws(GeneSymbol),                # 去除前后空格
    UniprotID = as.character(UniprotID),
    UniprotID = trimws(UniprotID)
  ) %>%
  # 移除重复的基因（如果有的话，保留FoldChange最高的）
  dplyr::group_by(GeneSymbol, Treatment) %>%
  dplyr::slice_max(FoldChange, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

cat("Data cleaning complete!\n")
cat("Before cleaning:", nrow(all_significant_proteins), "proteins\n")
cat("After cleaning:", nrow(all_significant_proteins_clean), "proteins\n")

# 检查清理结果
if (any(is.na(all_significant_proteins_clean$GeneSymbol))) {
  cat("Warning: Still found NA in GeneSymbol after cleaning\n")
  na_rows <- which(is.na(all_significant_proteins_clean$GeneSymbol))
  cat("NA rows:", paste(na_rows, collapse = ", "), "\n")
}

# 继续进行注释
annotated_proteins <- all_significant_proteins_clean %>%
  dplyr::mutate(
    # GO annotation matching - 确保匹配结果是逻辑值
    ER_match = GeneSymbol %in% df_er_prots_all$SYMBOL,
    Golgi_match = GeneSymbol %in% df_golgi_prots_all$SYMBOL,
    Mito_match = GeneSymbol %in% df_mito_prots_all$SYMBOL,
    Autophagy_match = GeneSymbol %in% df_autophagy_prots_all$SYMBOL,
    # LIR annotation matching
    LIR_match = GeneSymbol %in% df_lir_proteins$SYMBOL,
    # Any GO match
    Any_GO_match = ER_match | Golgi_match | Mito_match | Autophagy_match
  ) %>%
  # Add LIR detailed information
  dplyr::left_join(
    df_lir_proteins %>% dplyr::select(SYMBOL, LIR_Type, Confidence, anchor_no),
    by = c("GeneSymbol" = "SYMBOL")
  ) %>%
  dplyr::mutate(
    # 确保LIR相关字段没有问题
    LIR_Type = ifelse(is.na(LIR_Type), "No_LIR", as.character(LIR_Type)),
    Confidence = ifelse(is.na(Confidence), "None", as.character(Confidence)),
    anchor_no = ifelse(is.na(anchor_no), 0, as.numeric(anchor_no)),
    # 确保所有逻辑字段都是逻辑值
    ER_match = as.logical(ER_match),
    Golgi_match = as.logical(Golgi_match),
    Mito_match = as.logical(Mito_match),
    Autophagy_match = as.logical(Autophagy_match),
    LIR_match = as.logical(LIR_match),
    Any_GO_match = as.logical(Any_GO_match)
  )

# 最终检查
cat("Final annotation check:\n")
cat("- Total annotated proteins:", nrow(annotated_proteins), "\n")
cat("- Proteins with NA GeneSymbol:", sum(is.na(annotated_proteins$GeneSymbol)), "\n")
cat("- ER matches:", sum(annotated_proteins$ER_match, na.rm = TRUE), "\n")
cat("- Golgi matches:", sum(annotated_proteins$Golgi_match, na.rm = TRUE), "\n")
cat("- Mitochondrial matches:", sum(annotated_proteins$Mito_match, na.rm = TRUE), "\n")
cat("- Autophagy matches:", sum(annotated_proteins$Autophagy_match, na.rm = TRUE), "\n")
cat("- LIR matches:", sum(annotated_proteins$LIR_match, na.rm = TRUE), "\n")

# Display LIR matched proteins (排除NA)
lir_matches <- annotated_proteins %>% 
  dplyr::filter(LIR_match == TRUE, !is.na(GeneSymbol))

if (nrow(lir_matches) > 0) {
  cat("\nFound LIR proteins:\n")
  print(lir_matches %>% 
          dplyr::select(GeneSymbol, Treatment, FoldChange, LIR_Type, Confidence, anchor_no) %>%
          dplyr::arrange(desc(anchor_no)))
} else {
  cat("\nNo LIR proteins found\n")
}

# ============= Step 7: Create Statistics =============

# Statistics by treatment group
go_lir_summary_by_treatment <- annotated_proteins %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
    Total_Proteins = n(),
    ER_proteins = sum(ER_match),
    Golgi_proteins = sum(Golgi_match),
    Mito_proteins = sum(Mito_match),
    Autophagy_proteins = sum(Autophagy_match),
    LIR_proteins = sum(LIR_match),
    LIR_High_Confidence = sum(LIR_match & Confidence == "High"),
    LIR_Medium_Confidence = sum(LIR_match & Confidence == "Medium"),
    LIR_Low_Confidence = sum(LIR_match & Confidence == "Low"),
    No_Annotation = sum(!Any_GO_match & !LIR_match),
    .groups = 'drop'
  )

cat("\n=== GO and LIR Annotation Statistics by Treatment ===\n")
print(go_lir_summary_by_treatment)

# ============= Step 8: Create Compact Linear-Style Radial Plot Function =============

create_go_lir_radial_plot_en <- function(treatment_name, treatment_data) {
  
  # 首先清理数据 - 移除NA和空的基因符号
  treatment_data_clean <- treatment_data %>%
    dplyr::filter(
      !is.na(GeneSymbol),           # 移除NA
      GeneSymbol != "",             # 移除空字符串
      GeneSymbol != "NA",           # 移除字符串"NA"
      !is.na(FoldChange),           # 移除FoldChange为NA的行
      !is.infinite(FoldChange)      # 移除无限值
    ) %>%
    dplyr::mutate(
      GeneSymbol = as.character(GeneSymbol)  # 确保是字符类型
    )
  
  # 检查清理后是否还有数据
  if (nrow(treatment_data_clean) == 0) {
    cat(paste("Skipping", treatment_name, "- no valid proteins after data cleaning\n"))
    return(NULL)
  }
  
  # Sort by fold change in descending order
  treatment_data_ordered <- treatment_data_clean %>%
    dplyr::arrange(desc(FoldChange)) %>%
    dplyr::mutate(
      GeneSymbol = factor(GeneSymbol, levels = GeneSymbol),
      position = row_number()
    )
  
  cat(paste("Processing", treatment_name, "- valid proteins:", nrow(treatment_data_ordered), "\n"))
  
  # 检查是否有重复的基因符号
  duplicated_genes <- treatment_data_ordered$GeneSymbol[duplicated(treatment_data_ordered$GeneSymbol)]
  if (length(duplicated_genes) > 0) {
    cat("Warning: Found duplicated gene symbols:", paste(duplicated_genes, collapse = ", "), "\n")
    # 如果有重复，保留FoldChange最高的
    treatment_data_ordered <- treatment_data_ordered %>%
      dplyr::group_by(GeneSymbol) %>%
      dplyr::slice_max(FoldChange, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(FoldChange)) %>%
      dplyr::mutate(position = row_number())
  }
  
  # Prepare GO annotation plot data
  go_plot_data <- treatment_data_ordered %>%
    dplyr::select(GeneSymbol, position, ER_match, Golgi_match, Mito_match, Autophagy_match) %>%
    # 确保逻辑值正确
    dplyr::mutate(
      ER_match = as.logical(ER_match),
      Golgi_match = as.logical(Golgi_match),
      Mito_match = as.logical(Mito_match),
      Autophagy_match = as.logical(Autophagy_match)
    ) %>%
    # 替换NA为FALSE
    dplyr::mutate(
      ER_match = ifelse(is.na(ER_match), FALSE, ER_match),
      Golgi_match = ifelse(is.na(Golgi_match), FALSE, Golgi_match),
      Mito_match = ifelse(is.na(Mito_match), FALSE, Mito_match),
      Autophagy_match = ifelse(is.na(Autophagy_match), FALSE, Autophagy_match)
    ) %>%
    tidyr::pivot_longer(cols = c(ER_match, Golgi_match, Mito_match, Autophagy_match),
                        names_to = "GO_Category", values_to = "Match") %>%
    dplyr::filter(Match == TRUE) %>%  # 只保留TRUE的匹配
    dplyr::mutate(
      GO_Category = dplyr::case_when(
        GO_Category == "ER_match" ~ "ER",
        GO_Category == "Golgi_match" ~ "Golgi", 
        GO_Category == "Mito_match" ~ "Mitochondrion",
        GO_Category == "Autophagy_match" ~ "Autophagy"
      ),
      GO_Category = factor(GO_Category, levels = c("ER", "Golgi", "Mitochondrion", "Autophagy")),
      # 超紧凑的层布局，进一步减小层间距
      y_position = dplyr::case_when(
        GO_Category == "ER" ~ 0.6,
        GO_Category == "Golgi" ~ 1.0,
        GO_Category == "Mitochondrion" ~ 1.4,
        GO_Category == "Autophagy" ~ 1.8,
        TRUE ~ 0.6
      )
    ) %>%
    # 再次过滤掉任何可能的NA
    dplyr::filter(!is.na(GeneSymbol), !is.na(position), !is.na(GO_Category))
  
  # Prepare LIR plot data - outermost layer
  lir_plot_data <- treatment_data_ordered %>%
    dplyr::filter(
      LIR_match == TRUE,           # 只保留TRUE的LIR匹配
      !is.na(GeneSymbol),          # 确保基因符号不是NA
      !is.na(position)             # 确保位置不是NA
    ) %>%
    dplyr::mutate(
      y_position = 2.4,  # 更靠近其他层的LIR层位置
      # 确保Confidence和anchor_no不是NA
      Confidence = ifelse(is.na(Confidence), "None", as.character(Confidence)),
      anchor_no = ifelse(is.na(anchor_no), 0, as.numeric(anchor_no)),
      # Color based on confidence - 使用更类似第一张图的颜色
      LIR_color = dplyr::case_when(
        Confidence == "High" ~ "#008080",          # 深青色
        Confidence == "Medium" ~ "#20B2AA",        # 浅海绿
        Confidence == "Low" ~ "#48D1CC",           # 中等绿松石色
        Confidence == "Very_Low" ~ "#AFEEEE",      # 苍白绿松石色
        TRUE ~ "#CCCCCC"                           # 灰色
      ),
      # Size based on anchor_no - 更明显的大小差异
      LIR_size = dplyr::case_when(
        anchor_no >= 20 ~ 6,                       
        anchor_no >= 15 ~ 5,                       
        anchor_no >= 10 ~ 4,                       
        anchor_no >= 5 ~ 3,                        
        TRUE ~ 2                                   
      )
    )
  
  # Statistics
  n_total <- nrow(treatment_data_ordered)
  n_lir <- nrow(lir_plot_data)
  n_go <- length(unique(go_plot_data$position))
  
  # 调试信息
  cat(paste("  Total proteins:", n_total, "\n"))
  cat(paste("  GO matches:", n_go, "\n"))
  cat(paste("  LIR matches:", n_lir, "\n"))
  cat(paste("  Highest FC protein:", treatment_data_ordered$GeneSymbol[1], 
            "FC =", round(treatment_data_ordered$FoldChange[1], 2), "\n"))
  
  # Create basic radial plot - 更像第一张图的样式
  p <- ggplot() +
    # Add background ring - all proteins (更小更紧凑)
    geom_point(data = treatment_data_ordered,
               aes(x = position, y = 0), 
               color = "lightgray", size = 0.5, alpha = 0.3) +
    
    # Add GO annotation points (更小更紧凑的点)
    {
      if (nrow(go_plot_data) > 0) {
        geom_point(data = go_plot_data,
                   aes(x = position, y = y_position, color = GO_Category), 
                   size = 2, alpha = 0.9)
      }
    } +
    
    # Add LIR protein points (稍小的LIR点)
    {
      if (nrow(lir_plot_data) > 0) {
        geom_point(data = lir_plot_data,
                   aes(x = position, y = y_position, size = LIR_size), 
                   color = lir_plot_data$LIR_color, alpha = 0.9)
      }
    } +
    
    # Add LIR protein points (稍小的LIR点)
    {
      if (nrow(lir_plot_data) > 0) {
        geom_point(data = lir_plot_data,
                   aes(x = position, y = y_position, size = LIR_size), 
                   color = lir_plot_data$LIR_color, alpha = 0.9)
      }
    } +
    
    # Use polar coordinates - 创建有缺口的弧形 (240度弧，120度缺口)
    coord_polar(theta = "x", start = -pi/3, direction = 1) +
    
    # Set GO annotation colors - 更鲜明的对比色
    scale_color_manual(name = "GO Annotation",
                       values = c("ER" = "#E31A1C",           # 红色
                                  "Golgi" = "#1F78B4",        # 蓝色  
                                  "Mitochondrion" = "#33A02C", # 绿色
                                  "Autophagy" = "#FF7F00")) +  # 橙色
    
    # Set LIR point size
    scale_size_identity() +
    
    # Set y-axis - 拉近LIR层和蛋白质名称层的距离
    scale_y_continuous(
      breaks = c(0, 0.6, 1.0, 1.4, 1.8, 2.4), 
      labels = c("All", "1_ER", "2_Golgi", "3_Mito", "4_Autophagy", "5_LIR"),
      limits = c(-0.1, 2.7)  # 缩小上限，使整体更紧凑
    ) +
    
    # Set x-axis - 显示所有蛋白质名称，更大字体
    scale_x_continuous(
      breaks = 1:n_total,  # 显示所有蛋白质位置
      labels = as.character(treatment_data_ordered$GeneSymbol),  # 显示所有蛋白质名称
      limits = c(1, n_total)
    ) +
    
    # Theme settings - 超紧凑主题，更大字体
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "gray20", face = "bold"),  # 更大更粗的字体
      axis.text.y = element_text(size = 9, color = "gray30", face = "bold"),
      panel.grid.major = element_line(color = "gray88", size = 0.2),
      panel.grid.minor = element_line(color = "gray94", size = 0.1),
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, color = "gray20"),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 6, hjust = 0.5, color = "gray50"),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      plot.margin = margin(2, 2, 2, 2, "pt"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    
    # Titles and captions
    labs(
      title = paste("SL18", treatment_name, "- GO Annotation & LIR Analysis"),
      subtitle = paste("Significantly Upregulated Proteins (n =", n_total, ") - Ordered by Fold Change"),
      caption = paste("GO matches:", n_go, "| LIR matches:", n_lir, "| Ultra-compact layout | Ordered by fold change (highest→lowest)"),
      x = "",
      y = ""
    ) +
    
    # Add legend
    guides(
      color = guide_legend(
        title = "GO Annotation", 
        override.aes = list(size = 3),
        title.position = "top",
        title.hjust = 0.5
      ),
      size = "none"  
    )
  
  return(p)
}

# ============= Step 9: Generate All English Plots =============

# Create plots directory in current working directory
plots_dir <- "plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
  cat("Created plots directory at:", file.path(getwd(), plots_dir), "\n")
} else {
  cat("Using existing plots directory at:", file.path(getwd(), plots_dir), "\n")
}

cat("Current working directory:", getwd(), "\n")
cat("All plots will be saved in:", file.path(getwd(), plots_dir), "\n\n")

cat("Generating English plots...\n")

# 1. Save radial plots for each treatment group
cat("Saving radial plots...\n")
for (treat in unique(annotated_proteins$Treatment)) {
  treatment_data <- annotated_proteins %>%
    dplyr::filter(Treatment == treat)
  
  if (nrow(treatment_data) > 0) {
    cat("Processing", treat, "group...\n")
    
    p_radial <- create_go_lir_radial_plot_en(treat, treatment_data)
    
    if (!is.null(p_radial)) {
      filename <- file.path(plots_dir, paste0("SL18_", treat, "_GO_LIR_radial_plot_EN.png"))
      ggsave(filename, plot = p_radial, width = 14, height = 12, dpi = 300)
      cat("  ✓ Saved:", filename, "\n")
    }
  }
}

# 2. Save summary statistics plot
cat("Saving summary statistics plot...\n")
go_lir_summary_long <- go_lir_summary_by_treatment %>%
  tidyr::pivot_longer(cols = c(ER_proteins, Golgi_proteins, Mito_proteins, 
                               Autophagy_proteins, LIR_proteins, No_Annotation),
                      names_to = "Annotation_Category", values_to = "Count") %>%
  dplyr::mutate(Annotation_Category = dplyr::case_when(
    Annotation_Category == "ER_proteins" ~ "ER",
    Annotation_Category == "Golgi_proteins" ~ "Golgi",
    Annotation_Category == "Mito_proteins" ~ "Mitochondrion",
    Annotation_Category == "Autophagy_proteins" ~ "Autophagy",
    Annotation_Category == "LIR_proteins" ~ "LIR (ATG8 Interaction)",
    Annotation_Category == "No_Annotation" ~ "No Annotation"
  ))

p_summary <- go_lir_summary_long %>%
  ggplot(aes(x = Treatment, y = Count, fill = Annotation_Category)) +
  geom_col(position = "stack", alpha = 0.8) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")), 
            position = position_stack(vjust = 0.5), size = 3) +
  labs(
    title = "SL18 GO Annotation and LIR Distribution by Treatment",
    subtitle = "Significantly Upregulated Proteins - GO Annotation and ATG8 Interaction Analysis",
    x = "Treatment Group",
    y = "Number of Proteins",
    fill = "Annotation Category"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  scale_fill_manual(values = c(
    "ER" = "#E31A1C",
    "Golgi" = "#1F78B4", 
    "Mitochondrion" = "#33A02C",
    "Autophagy" = "#FF7F00",
    "LIR (ATG8 Interaction)" = "#6A3D9A",
    "No Annotation" = "#CCCCCC"
  ))

# Save summary plot
filename_summary <- file.path(plots_dir, "SL18_GO_LIR_summary_plot_EN.png")
ggsave(filename_summary, plot = p_summary, width = 12, height = 8, dpi = 300)
cat("  ✓ Saved:", filename_summary, "\n")

# 3. Display final summary
cat("\n=== Plot Generation Complete ===\n")
cat("All plots saved in:", file.path(getwd(), plots_dir), "\n")
cat("Generated files:\n")

# List all generated files
plot_files <- list.files(plots_dir, pattern = "\\.png$", full.names = FALSE)
for (i in seq_along(plot_files)) {
  cat(sprintf("  %d. %s\n", i, plot_files[i]))
}

cat("\nTotal plots generated:", length(plot_files), "\n")

# Optional: Show file sizes
cat("\nFile sizes:\n")
for (file in plot_files) {
  full_path <- file.path(plots_dir, file)
  if (file.exists(full_path)) {
    size_mb <- round(file.info(full_path)$size / 1024 / 1024, 2)
    cat(sprintf("  %s: %.2f MB\n", file, size_mb))
  }
}