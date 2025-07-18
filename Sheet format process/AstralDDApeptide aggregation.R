# 蛋白质数据合并处理脚本,用于将AstralDDA的肽段数据合并来增强信号，之后进行过滤，并按照符合TMTviewer输入格式的方式进行输出
# 按ProteinId合并行，对定量数据求和，添加肽段计数

library(dplyr)
library(readr)

# 定义处理函数
process_protein_data <- function(input_file, output_file) {
  
  # 读取TSV文件
  cat("正在读取文件:", input_file, "\n")
  data <- read_tsv(input_file, show_col_types = FALSE)
  
  # 检查必需的列是否存在
  required_cols <- c("ProteinId", "GeneSymbol", "Description", "GroupId")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("缺少必需的列: ", paste(missing_cols, collapse = ", "))
  }
  
  # 检查Collapsed列是否存在，如果不存在则创建默认值
  if (!"Collapsed" %in% colnames(data)) {
    cat("警告: Collapsed列不存在，将创建默认值\n")
    data$Collapsed <- "C"  # 设置默认值
  }
  
  # 检查Unique/Razor列是否存在
  if (!"Unique/Razor" %in% colnames(data)) {
    cat("警告: Unique/Razor列不存在，Parsimony将基于肽段数量设置\n")
    use_ur_info <- FALSE
  } else {
    use_ur_info <- TRUE
    cat("找到Unique/Razor列，将用于设置Parsimony信息\n")
  }
  
  # 定义需要合并求和的rq列
  rq_columns <- c("rq_126_sn", "rq_127n_sn", "rq_127c_sn", "rq_128n_sn", 
                  "rq_128c_sn", "rq_129n_sn", "rq_129c_sn", "rq_130n_sn", 
                  "rq_130c_sn", "rq_131_sn", "rq_131c_sn", "rq_132n_sn", 
                  "rq_132c_sn", "rq_133n_sn", "rq_133c_sn", "rq_134n_sn", 
                  "rq_134c_sn", "rq_135n_sn")
  
  # 检查rq列是否存在
  existing_rq_cols <- intersect(rq_columns, colnames(data))
  missing_rq_cols <- setdiff(rq_columns, colnames(data))
  
  if (length(missing_rq_cols) > 0) {
    cat("警告: 以下rq列不存在，将被跳过:", paste(missing_rq_cols, collapse = ", "), "\n")
  }
  
  cat("找到", length(existing_rq_cols), "个rq列用于合并\n")
  
  # 按ProteinId分组并合并数据
  cat("正在按ProteinId分组合并数据...\n")
  
  if (use_ur_info) {
    # 如果有Unique/Razor信息，根据该信息设置Parsimony
    merged_data <- data %>%
      group_by(ProteinId, GeneSymbol, Description, GroupId) %>%
      summarise(
        # 计算每个蛋白质的肽段数量
        `Number of peptides` = n(),
        
        # 保留Collapsed列（取第一个值，通常相同蛋白质的该值应该一致）
        Collapsed = first(Collapsed, na_rm = TRUE),
        
        # 根据Unique/Razor信息设置Parsimony
        Parsimony = case_when(
          # 如果所有肽段都是Unique
          all(`Unique/Razor` == "U", na.rm = TRUE) ~ "Unique",
          # 如果包含Unique肽段
          any(`Unique/Razor` == "U", na.rm = TRUE) ~ "Mixed",
          # 如果所有肽段都是Razor
          all(`Unique/Razor` == "R", na.rm = TRUE) ~ "Razor",
          # 如果包含Shared肽段
          any(`Unique/Razor` == "S", na.rm = TRUE) ~ "Shared",
          # 默认情况
          TRUE ~ "Unknown"
        ),
        
        # 对所有存在的rq列求和，处理NA值
        across(all_of(existing_rq_cols), ~ sum(.x, na.rm = TRUE)),
        
        .groups = "drop"
      ) %>%
      # 按照要求的顺序重新排列列
      select(ProteinId, GeneSymbol, Description, GroupId, `Number of peptides`, 
             Collapsed, Parsimony, all_of(existing_rq_cols))
  } else {
    # 如果没有Unique/Razor信息，使用原来的逻辑
    merged_data <- data %>%
      group_by(ProteinId, GeneSymbol, Description, GroupId) %>%
      summarise(
        # 计算每个蛋白质的肽段数量
        `Number of peptides` = n(),
        
        # 保留Collapsed列（取第一个值，通常相同蛋白质的该值应该一致）
        Collapsed = first(Collapsed, na_rm = TRUE),
        
        # 添加Parsimony列（设置为合并后的标识，或者可以根据需要修改逻辑）
        Parsimony = if_else(n() > 1, "Merged", "Single"),
        
        # 对所有存在的rq列求和，处理NA值
        across(all_of(existing_rq_cols), ~ sum(.x, na.rm = TRUE)),
        
        .groups = "drop"
      ) %>%
      # 按照要求的顺序重新排列列
      select(ProteinId, GeneSymbol, Description, GroupId, `Number of peptides`, 
             Collapsed, Parsimony, all_of(existing_rq_cols))
  }
  
  # 显示处理结果统计
  cat("处理完成!\n")
  cat("原始数据行数:", nrow(data), "\n")
  cat("合并后行数:", nrow(merged_data), "\n")
  cat("唯一蛋白质数量:", length(unique(data$ProteinId)), "\n")
  
  # 如果使用了Unique/Razor信息，显示Parsimony分布
  if (use_ur_info) {
    cat("\nParsimony分布:\n")
    parsimony_counts <- table(merged_data$Parsimony)
    for (i in 1:length(parsimony_counts)) {
      cat("  ", names(parsimony_counts)[i], ":", parsimony_counts[i], "\n")
    }
  }
  
  # 显示每个蛋白质的肽段数量分布
  peptide_counts <- merged_data$`Number of peptides`
  cat("\n肽段数量统计:\n")
  cat("  最小值:", min(peptide_counts), "\n")
  cat("  最大值:", max(peptide_counts), "\n")
  cat("  平均值:", round(mean(peptide_counts), 2), "\n")
  cat("  中位数:", median(peptide_counts), "\n")
  
  # 写入第一个输出文件（完整数据）
  cat("正在写入完整数据输出文件:", output_file, "\n")
  
  # 重命名列标题，添加空格
  output_data <- merged_data %>%
    rename(
      `Protein Id` = ProteinId,
      `Gene Symbol` = GeneSymbol,
      `Group ID` = GroupId
    )
  
  # 创建scaled列 - 每个rq列除以所有rq列的总和再乘以100
  cat("正在计算scaled列...\n")
  
  # 计算每行所有rq列的总和
  output_data <- output_data %>%
    rowwise() %>%
    mutate(
      total_rq_sum = sum(c_across(all_of(existing_rq_cols)), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # 创建scaled列，使用符合正则表达式的列名格式
  for (col in existing_rq_cols) {
    # 将 rq_126_sn 转换为 126_sn_scaled 格式
    new_col_name <- gsub("rq_", "", col)  # 移除 "rq_" 前缀
    scaled_col_name <- paste0(new_col_name, "_scaled")
    
    output_data[[scaled_col_name]] <- ifelse(
      output_data$total_rq_sum == 0, 
      0, 
      (output_data[[col]] / output_data$total_rq_sum) * 100
    )
  }
  
  # 移除临时的总和列
  output_data <- output_data %>%
    select(-total_rq_sum)
  
  # 重命名sum列，使用符合正则表达式的列名格式
  rq_rename_list <- c()
  for (col in existing_rq_cols) {
    new_col_name <- gsub("rq_", "", col)  # 移除 "rq_" 前缀
    sum_col_name <- paste0(new_col_name, "_sum")
    rq_rename_list[sum_col_name] <- col
  }
  
  output_data <- output_data %>%
    rename(all_of(rq_rename_list))
  
  # 重新排列列的顺序：基本信息 + scaled列 + sum列
  scaled_cols <- paste0(gsub("rq_", "", existing_rq_cols), "_scaled")
  sum_cols <- paste0(gsub("rq_", "", existing_rq_cols), "_sum")
  
  output_data <- output_data %>%
    select(`Protein Id`, `Gene Symbol`, Description, `Group ID`, `Number of peptides`, 
           Collapsed, Parsimony, all_of(scaled_cols), all_of(sum_cols))
  
  write_tsv(output_data, output_file)
  
  cat("完整数据文件处理完成！输出文件已保存为:", output_file, "\n")
  
  # 创建过滤后的数据集
  cat("\n开始进行数据过滤...\n")
  
  # 计算每行所有rq列的总和
  merged_data_with_sum <- merged_data %>%
    rowwise() %>%
    mutate(
      rq_total_sum = sum(c_across(all_of(existing_rq_cols)), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # 过滤掉总和小于1400的行
  filtered_data <- merged_data_with_sum %>%
    filter(rq_total_sum >= 1400) %>%
    select(-rq_total_sum)  # 移除临时的总和列
  
  # 计算过滤统计
  original_count <- nrow(merged_data)
  filtered_count <- nrow(filtered_data)
  removed_count <- original_count - filtered_count
  removed_percentage <- round((removed_count / original_count) * 100, 2)
  
  # 显示过滤统计
  cat("数据过滤统计:\n")
  cat("  原始蛋白质数量:", original_count, "\n")
  cat("  过滤后蛋白质数量:", filtered_count, "\n")
  cat("  被移除的蛋白质数量:", removed_count, "\n")
  cat("  被移除的百分比:", removed_percentage, "%\n")
  
  # 生成过滤后的文件名
  filtered_output_file <- gsub("\\.tsv$", "_filtered.tsv", output_file)
  
  # 写入过滤后的数据文件
  cat("正在写入过滤后数据文件:", filtered_output_file, "\n")
  
  # 重命名列标题，添加空格
  filtered_output_data <- filtered_data %>%
    rename(
      `Protein Id` = ProteinId,
      `Gene Symbol` = GeneSymbol,
      `Group ID` = GroupId
    )
  
  # 创建scaled列 - 每个rq列除以所有rq列的总和再乘以100
  cat("正在计算过滤数据的scaled列...\n")
  
  # 计算每行所有rq列的总和
  filtered_output_data <- filtered_output_data %>%
    rowwise() %>%
    mutate(
      total_rq_sum = sum(c_across(all_of(existing_rq_cols)), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # 创建scaled列，使用符合正则表达式的列名格式
  for (col in existing_rq_cols) {
    # 将 rq_126_sn 转换为 126_sn_scaled 格式
    new_col_name <- gsub("rq_", "", col)  # 移除 "rq_" 前缀
    scaled_col_name <- paste0(new_col_name, "_scaled")
    
    filtered_output_data[[scaled_col_name]] <- ifelse(
      filtered_output_data$total_rq_sum == 0, 
      0, 
      (filtered_output_data[[col]] / filtered_output_data$total_rq_sum) * 100
    )
  }
  
  # 移除临时的总和列
  filtered_output_data <- filtered_output_data %>%
    select(-total_rq_sum)
  
  # 重命名sum列，使用符合正则表达式的列名格式
  rq_rename_list <- c()
  for (col in existing_rq_cols) {
    new_col_name <- gsub("rq_", "", col)  # 移除 "rq_" 前缀
    sum_col_name <- paste0(new_col_name, "_sum")
    rq_rename_list[sum_col_name] <- col
  }
  
  filtered_output_data <- filtered_output_data %>%
    rename(all_of(rq_rename_list))
  
  # 重新排列列的顺序：基本信息 + scaled列 + sum列
  scaled_cols <- paste0(gsub("rq_", "", existing_rq_cols), "_scaled")
  sum_cols <- paste0(gsub("rq_", "", existing_rq_cols), "_sum")
  
  filtered_output_data <- filtered_output_data %>%
    select(`Protein Id`, `Gene Symbol`, Description, `Group ID`, `Number of peptides`, 
           Collapsed, Parsimony, all_of(scaled_cols), all_of(sum_cols))
  
  write_tsv(filtered_output_data, filtered_output_file)
  
  cat("过滤后数据文件已保存为:", filtered_output_file, "\n")
  
  # 创建并保存处理日志
  log_info <- data.frame(
    处理时间 = as.character(Sys.time()),
    输入文件 = input_file,
    完整输出文件 = output_file,
    过滤输出文件 = filtered_output_file,
    原始行数 = nrow(data),
    合并后行数 = nrow(merged_data),
    过滤后行数 = nrow(filtered_data),
    唯一蛋白质数 = nrow(merged_data),
    过滤后蛋白质数 = nrow(filtered_data),
    被移除蛋白质数 = removed_count,
    被移除百分比 = paste0(removed_percentage, "%")
  )
  
  write_tsv(log_info, "processing_log.tsv")
  cat("\n处理日志已保存为: processing_log.tsv\n")
  
  # 创建详细的过滤报告
  filter_report <- data.frame(
    项目 = c("处理时间", "输入文件", "完整输出文件", "过滤输出文件", 
           "原始数据行数", "合并后蛋白质数量", "过滤后蛋白质数量", 
           "被移除蛋白质数量", "被移除百分比", "过滤标准"),
    值 = c(as.character(Sys.time()), input_file, output_file, filtered_output_file,
          nrow(data), nrow(merged_data), nrow(filtered_data),
          removed_count, paste0(removed_percentage, "%"), 
          "rq列总和 >= 1400")
  )
  
  write_tsv(filter_report, "filter_report.tsv")
  cat("过滤详细报告已保存为: filter_report.tsv\n")
  
  # 返回包含两个数据集的列表
  return(list(
    complete_data = merged_data,
    filtered_data = filtered_data,
    filter_stats = list(
      original_count = original_count,
      filtered_count = filtered_count,
      removed_count = removed_count,
      removed_percentage = removed_percentage
    )
  ))
}

# 使用示例
# 请根据您的实际文件路径修改以下路径
input_file <- "input_data.tsv"  # 输入TSV文件路径
output_file <- "merged_protein_data.tsv"  # 输出文件路径

# 执行处理
result <- process_protein_data(input_file, output_file)

# 显示前几行完整结果
cat("\n前5行完整数据预览:\n")
print(head(result$complete_data, 5))

# 显示前几行过滤后结果
cat("\n前5行过滤后数据预览:\n")
if (nrow(result$filtered_data) > 0) {
  print(head(result$filtered_data, 5))
} else {
  cat("没有数据通过过滤条件\n")
}

# 显示过滤统计摘要
cat("\n=== 数据处理摘要 ===\n")
cat("完整数据蛋白质数量:", result$filter_stats$original_count, "\n")
cat("过滤后蛋白质数量:", result$filter_stats$filtered_count, "\n")
cat("被移除蛋白质数量:", result$filter_stats$removed_count, "\n")
cat("被移除百分比:", result$filter_stats$removed_percentage, "%\n")

# 可选：显示特定蛋白质的详细信息
if (nrow(result$complete_data) > 0) {
  cat("\n示例蛋白质详细信息（完整数据）:\n")
  sample_protein <- result$complete_data[1, ]
  print(sample_protein)
} else {
  cat("\n没有找到处理结果\n")
}