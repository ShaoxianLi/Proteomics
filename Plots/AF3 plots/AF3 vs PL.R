# 清理环境并重新加载包
rm(list = ls())

# 加载必要的包（按正确顺序）
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(gridExtra)
library(ggrepel)  # 用于文本标签

# 解决select函数冲突
select <- dplyr::select

# 设置文件路径
interaction_file <- "C:/Users/argo_/OneDrive - UMass Chan Medical School/Desktop/Result/plot/SL18_20250625_APEX_LC3B diff autophagy treatment on Eclipse/LC3B_interaction_summary.csv"
proteomics_file <- "C:/Users/argo_/OneDrive - UMass Chan Medical School/Desktop/Result/plot/SL18_20250625_APEX_LC3B diff autophagy treatment on Eclipse/SL18_Live_pq_153 modified.csv"

# 读取数据
interaction_data <- read_csv(interaction_file, show_col_types = FALSE)
proteomics_data <- read_csv(proteomics_file, show_col_types = FALSE)

print("=== 数据读取完成 ===")
print(paste("Interaction data rows:", nrow(interaction_data)))
print(paste("Proteomics data rows:", nrow(proteomics_data)))

# 检查列结构
print("Interaction data columns:")
print(colnames(interaction_data))
print("Proteomics data columns:")
print(colnames(proteomics_data))

# 检查前几行数据格式
print("Interaction data sample (first 3 rows, first 3 columns):")
print(head(interaction_data[, 1:min(3, ncol(interaction_data))], 3))
print("Proteomics data sample (first 3 rows, first 3 columns):")
print(head(proteomics_data[, 1:min(3, ncol(proteomics_data))], 3))

# =============================================================================
# UniProt ID 提取和匹配
# =============================================================================

# 从interaction数据的第二列提取UniProt ID
if (ncol(interaction_data) >= 2) {
  interaction_data$UniProt_ID_Int <- str_extract(interaction_data[[2]], "[POQ][0-9][A-Z0-9]{3}[0-9]")
  
  if (sum(!is.na(interaction_data$UniProt_ID_Int)) == 0) {
    interaction_data$UniProt_ID_Int <- str_extract(interaction_data[[1]], "[POQ][0-9][A-Z0-9]{3}[0-9]")
  }
} else {
  interaction_data$UniProt_ID_Int <- NA
}

# 从proteomics数据的第二列提取UniProt ID
if (ncol(proteomics_data) >= 2) {
  proteomics_data$UniProt_ID_Prot <- str_extract(proteomics_data[[2]], "[POQ][0-9][A-Z0-9]{3}[0-9]")
  
  if (sum(!is.na(proteomics_data$UniProt_ID_Prot)) == 0) {
    proteomics_data$UniProt_ID_Prot <- str_extract(proteomics_data[[1]], "[POQ][0-9][A-Z0-9]{3}[0-9]")
  }
} else {
  proteomics_data$UniProt_ID_Prot <- NA
}

# 备用：提取基因符号
interaction_data$Gene_Symbol_Int <- str_remove(interaction_data[[1]], "_HUMAN$")
proteomics_data$Gene_Symbol_Prot <- str_extract(proteomics_data[[1]], "^[^_]+")

print("=== UniProt ID 提取完成 ===")

# =============================================================================
# 数据匹配
# =============================================================================

# 通过UniProt ID匹配
interaction_with_uniprot <- interaction_data %>%
  filter(!is.na(UniProt_ID_Int))

proteomics_with_uniprot <- proteomics_data %>%
  filter(!is.na(UniProt_ID_Prot))

uniprot_matches <- interaction_with_uniprot %>%
  inner_join(proteomics_with_uniprot, 
             by = c("UniProt_ID_Int" = "UniProt_ID_Prot")) %>%
  mutate(Match_Type = "UniProt_ID")

# 对于未匹配的，尝试基因符号匹配
if (nrow(uniprot_matches) > 0) {
  matched_interaction_names <- uniprot_matches[[colnames(interaction_data)[1]]]
} else {
  matched_interaction_names <- character(0)
}

remaining_interaction <- interaction_data %>%
  filter(!(!!sym(colnames(interaction_data)[1]) %in% matched_interaction_names))

gene_matches <- remaining_interaction %>%
  inner_join(proteomics_data, 
             by = c("Gene_Symbol_Int" = "Gene_Symbol_Prot")) %>%
  mutate(Match_Type = "Gene_Symbol")

# 合并匹配结果
if (nrow(uniprot_matches) > 0 && nrow(gene_matches) > 0) {
  common_cols <- intersect(colnames(uniprot_matches), colnames(gene_matches))
  matched_data <- bind_rows(
    uniprot_matches %>% select(all_of(common_cols)),
    gene_matches %>% select(all_of(common_cols))
  )
} else if (nrow(uniprot_matches) > 0) {
  matched_data <- uniprot_matches
} else if (nrow(gene_matches) > 0) {
  matched_data <- gene_matches
} else {
  stop("没有匹配到任何蛋白质，请检查数据格式")
}

print("=== 匹配结果 ===")
print(paste("总匹配数:", nrow(matched_data)))
print(paste("匹配率:", round(nrow(matched_data)/nrow(interaction_data)*100, 1), "%"))

# =============================================================================
# 数值数据处理和统计分析
# =============================================================================

# 检查数值列
treatments <- c("Torin1_Baf", "MLSA5_Baf", "OA_Baf", "Mon_Baf", "Tunica_Baf")
control_cols <- c("Ctr_1", "Ctr_2", "Ctr_3")

# 确保数值列是数值类型
numeric_cols <- c(control_cols, paste0(rep(treatments, each = 3), "_", rep(1:3, length(treatments))))
existing_cols <- intersect(numeric_cols, colnames(matched_data))

print("=== 检查数值列 ===")
print("Available numeric columns:")
print(existing_cols)

# 检查并修复数据类型问题
for (col in existing_cols) {
  col_class <- class(matched_data[[col]])
  if ("list" %in% col_class) {
    print(paste("Converting", col, "from list to vector"))
    matched_data[[col]] <- unlist(matched_data[[col]])
  }
  matched_data[[col]] <- as.numeric(matched_data[[col]])
}

# =============================================================================
# 为每个treatment分别计算统计结果
# =============================================================================

# 存储所有treatment的结果
all_treatment_results <- list()

for (treatment in treatments) {
  treatment_cols <- paste0(treatment, "_", 1:3)
  
  # 检查列是否存在
  available_treatment_cols <- intersect(treatment_cols, colnames(matched_data))
  available_control_cols <- intersect(control_cols, colnames(matched_data))
  
  print(paste("Processing treatment:", treatment))
  
  if (length(available_treatment_cols) >= 2 && length(available_control_cols) >= 2) {
    
    # 计算平均值和fold change
    temp_data <- matched_data
    
    # 计算control平均值
    control_matrix <- as.matrix(temp_data[, available_control_cols])
    temp_data$Control_Mean <- rowMeans(control_matrix, na.rm = TRUE)
    
    # 计算treatment平均值
    treatment_matrix <- as.matrix(temp_data[, available_treatment_cols])
    temp_data$Treatment_Mean <- rowMeans(treatment_matrix, na.rm = TRUE)
    
    # 计算fold change
    temp_stats <- temp_data %>%
      filter(Control_Mean > 0, Treatment_Mean > 0, 
             is.finite(Control_Mean), is.finite(Treatment_Mean)) %>%
      mutate(
        Fold_Change = Treatment_Mean / Control_Mean,
        Log2_Fold_Change = log2(Treatment_Mean / Control_Mean),
        Treatment_Group = treatment
      ) %>%
      filter(is.finite(Fold_Change), is.finite(Log2_Fold_Change))
    
    # 进行t检验
    temp_stats$P_Value <- NA
    
    for (i in 1:nrow(temp_stats)) {
      control_vals <- as.numeric(temp_stats[i, available_control_cols])
      treatment_vals <- as.numeric(temp_stats[i, available_treatment_cols])
      
      control_vals <- control_vals[!is.na(control_vals) & control_vals > 0]
      treatment_vals <- treatment_vals[!is.na(treatment_vals) & treatment_vals > 0]
      
      if (length(control_vals) >= 2 && length(treatment_vals) >= 2 && 
          var(control_vals) > 0 && var(treatment_vals) > 0) {
        tryCatch({
          t_test <- t.test(treatment_vals, control_vals, var.equal = TRUE)
          temp_stats$P_Value[i] <- t_test$p.value
        }, error = function(e) {
          temp_stats$P_Value[i] <- NA
        })
      }
    }
    
    # 定义显著性（使用p value < 0.05 和 log2FC > 0.32）
    temp_stats <- temp_stats %>%
      mutate(
        Significant = !is.na(P_Value) & P_Value < 0.05 & abs(Log2_Fold_Change) > 0.32,
        Direction = case_when(
          !is.na(P_Value) & P_Value < 0.05 & Log2_Fold_Change > 0.32 ~ "Upregulated",
          !is.na(P_Value) & P_Value < 0.05 & Log2_Fold_Change < -0.32 ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      )
    
    # 添加蛋白质索引（按Total_High_Conf_Contacts排序）
    temp_stats <- temp_stats %>%
      arrange(desc(Total_High_Conf_Contacts)) %>%
      mutate(Protein_Index = row_number())
    
    all_treatment_results[[treatment]] <- temp_stats
  }
}

print("=== 所有treatment统计分析完成 ===")

# =============================================================================
# 创建6个图表
# =============================================================================

# 定义颜色
treatment_colors <- c(
  "Upregulated" = "#FF6B6B",     # 红色 - 上调
  "Downregulated" = "#4169E1",   # 蓝色 - 下调
  "Not Significant" = "#95A5A6"  # 灰色 - 不显著
)

# 创建单个treatment的图表函数
create_treatment_plot <- function(data, treatment_name) {
  # 添加绘图顺序和蛋白质名称
  plot_data <- data %>%
    mutate(
      Plot_Order = case_when(
        Direction == "Not Significant" ~ 1,
        Direction == "Downregulated" ~ 2,
        Direction == "Upregulated" ~ 3
      ),
      # 提取蛋白质名称用于标签
      Protein_Name = if("Gene_Symbol_Int" %in% colnames(data)) {
        Gene_Symbol_Int
      } else {
        str_extract(colnames(data)[1], "^[^_]+")
      },
      # 只为显著的蛋白质添加标签
      Label = ifelse(Significant, Protein_Name, "")
    ) %>%
    arrange(Plot_Order, Protein_Index)
  
  ggplot(plot_data, aes(x = Protein_Index, y = Total_High_Conf_Contacts)) +
    geom_point(aes(color = Direction), size = 2.5, alpha = 0.8) +
    geom_text_repel(aes(label = Label), 
                    size = 2.5, 
                    max.overlaps = 20,
                    box.padding = 0.3,
                    point.padding = 0.3,
                    segment.size = 0.2,
                    show.legend = FALSE) +
    scale_color_manual(values = treatment_colors, 
                       name = "Expression Status") +
    labs(
      title = paste("LC3B Interactions:", treatment_name, "vs Control"),
      subtitle = paste("Significant: |Log2FC| > 0.32 & p < 0.05"),
      x = "Protein Rank (by Total High-Confidence Contacts)",
      y = "Total High-Confidence Contacts"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
}

# 创建前5个图（每个treatment vs control）
plot_list <- list()
for (i in 1:length(treatments)) {
  treatment <- treatments[i]
  if (treatment %in% names(all_treatment_results)) {
    plot_list[[i]] <- create_treatment_plot(all_treatment_results[[treatment]], treatment)
  }
}

# =============================================================================
# 创建汇总图（第6个图）
# =============================================================================

# 合并所有treatment结果，为每个蛋白质找到最显著的结果
# 首先确定用于分组的列
group_col <- if ("UniProt_ID_Int" %in% colnames(all_treatment_results[[1]])) {
  "UniProt_ID_Int"
} else if ("Gene_Symbol_Int" %in% colnames(all_treatment_results[[1]])) {
  "Gene_Symbol_Int"
} else {
  colnames(all_treatment_results[[1]])[1]
}

# 合并所有结果
all_results_combined <- do.call(rbind, all_treatment_results)

# 为汇总图准备数据：如果蛋白质在任何一个treatment中显著，就标记为显著
summary_data <- all_results_combined %>%
  group_by(!!sym(group_col)) %>%
  summarise(
    # 添加蛋白质名称
    Protein_Name = if("Gene_Symbol_Int" %in% colnames(all_results_combined)) {
      Gene_Symbol_Int[1]
    } else {
      colnames(all_results_combined)[1][1]
    },
    Total_High_Conf_Contacts = Total_High_Conf_Contacts[1],
    # 如果任何一个treatment中显著，就标记为显著
    Any_Significant = any(Significant, na.rm = TRUE),
    # 记录最显著的结果
    Best_P_Value = min(P_Value, na.rm = TRUE),
    Best_Treatment = Treatment_Group[which.min(P_Value)][1],
    Best_Log2FC = Log2_Fold_Change[which.min(P_Value)][1],
    Best_Direction = Direction[which.min(P_Value)][1],
    # 统计显著的treatment数量
    Significant_Treatments = sum(Significant, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(Total_High_Conf_Contacts)) %>%
  mutate(
    Protein_Index = row_number(),
    Summary_Status = case_when(
      Any_Significant & Best_Direction == "Upregulated" ~ "Upregulated in Any Treatment",
      Any_Significant & Best_Direction == "Downregulated" ~ "Downregulated in Any Treatment",
      TRUE ~ "Not Significant"
    ),
    # 创建绘图顺序：不显著的在最下层
    Plot_Order = case_when(
      Summary_Status == "Not Significant" ~ 1,
      Summary_Status == "Downregulated in Any Treatment" ~ 2,
      Summary_Status == "Upregulated in Any Treatment" ~ 3
    )
  )

# 定义汇总图的颜色
summary_colors <- c(
  "Upregulated in Any Treatment" = "#FF6B6B",     # 红色
  "Downregulated in Any Treatment" = "#4169E1",   # 蓝色
  "Not Significant" = "#95A5A6"                   # 灰色
)

# 创建汇总图
summary_plot_data <- summary_data %>%
  arrange(Plot_Order, Protein_Index) %>%
  mutate(
    # 只为显著的蛋白质添加标签
    Label = ifelse(Any_Significant, Protein_Name, "")
  )

summary_plot <- ggplot(summary_plot_data, aes(x = Protein_Index, y = Total_High_Conf_Contacts)) +
  geom_point(aes(color = Summary_Status), size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = Label), 
                  size = 2.5, 
                  max.overlaps = 25,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.size = 0.2,
                  show.legend = FALSE) +
  scale_color_manual(values = summary_colors, 
                     name = "Expression Status") +
  labs(
    title = "LC3B Interactions: Summary of All Treatments",
    subtitle = "Proteins significant in any treatment (|Log2FC| > 0.32 & p < 0.05)",
    x = "Protein Rank (by Total High-Confidence Contacts)",
    y = "Total High-Confidence Contacts",
    caption = paste("Total proteins analyzed:", nrow(summary_data),
                    "| Significant in any treatment:", sum(summary_data$Any_Significant))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 8)
  )

plot_list[[6]] <- summary_plot

# =============================================================================
# 显示和保存图表
# =============================================================================

# 显示所有图表
for (i in 1:length(plot_list)) {
  if (!is.null(plot_list[[i]])) {
    print(plot_list[[i]])
  }
}

# 保存单独的图表
for (i in 1:5) {
  if (!is.null(plot_list[[i]])) {
    filename <- paste0("LC3B_", treatments[i], "_vs_Control.png")
    ggsave(filename, plot = plot_list[[i]], width = 10, height = 6, dpi = 300)
    cat("Saved:", filename, "\n")
  }
}

# 保存汇总图
ggsave("LC3B_Summary_All_Treatments.png", plot = summary_plot, width = 10, height = 6, dpi = 300)
cat("Saved: LC3B_Summary_All_Treatments.png\n")

# 保存组合图（可选）
if (length(plot_list) >= 6) {
  # 创建2x3布局的组合图
  combined_plot <- grid.arrange(grobs = plot_list, ncol = 2, nrow = 3)
  ggsave("LC3B_All_Plots_Combined.png", combined_plot, width = 16, height = 18, dpi = 300)
  cat("Saved: LC3B_All_Plots_Combined.png\n")
}

# =============================================================================
# 输出统计摘要
# =============================================================================

cat("\n=== Analysis Summary ===\n")
cat("Total proteins matched:", nrow(matched_data), "\n\n")

# 每个treatment的统计
for (treatment in treatments) {
  if (treatment %in% names(all_treatment_results)) {
    data <- all_treatment_results[[treatment]]
    sig_up <- sum(data$Direction == "Upregulated", na.rm = TRUE)
    sig_down <- sum(data$Direction == "Downregulated", na.rm = TRUE)
    cat(sprintf("%s vs Control:\n", treatment))
    cat(sprintf("  Upregulated: %d proteins\n", sig_up))
    cat(sprintf("  Downregulated: %d proteins\n", sig_down))
    cat(sprintf("  Total significant: %d proteins\n\n", sig_up + sig_down))
  }
}

# 汇总统计
cat("Summary across all treatments:\n")
cat(sprintf("  Proteins significant in any treatment: %d\n", sum(summary_data$Any_Significant)))
cat(sprintf("  Upregulated in any treatment: %d\n", 
            sum(summary_data$Summary_Status == "Upregulated in Any Treatment")))
cat(sprintf("  Downregulated in any treatment: %d\n", 
            sum(summary_data$Summary_Status == "Downregulated in Any Treatment")))

# 保存结果数据
write_csv(summary_data, "LC3B_Summary_Results.csv")
write_csv(all_results_combined, "LC3B_All_Treatment_Results.csv")

cat("\nData files saved:\n")
cat("- LC3B_Summary_Results.csv\n")
cat("- LC3B_All_Treatment_Results.csv\n")

# 输出前10个最高相互作用的蛋白质及其在各treatment中的表现
top_proteins <- summary_data %>%
  slice_head(n = 10) %>%
  select(!!sym(group_col), Total_High_Conf_Contacts, Any_Significant, 
         Best_Treatment, Best_Direction, Significant_Treatments)

cat("\nTop 10 LC3B-interacting proteins:\n")
print(top_proteins)