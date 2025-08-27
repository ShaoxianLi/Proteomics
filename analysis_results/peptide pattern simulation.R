# Peptide Isotope Distribution Analyzer - Clean Version
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(dplyr)

# 氨基酸分子式数据库 (C, H, N, O, S)
amino_acids <- list(
  A = c(3, 7, 1, 2, 0), R = c(6, 14, 4, 2, 0), N = c(4, 8, 2, 3, 0),
  D = c(4, 7, 1, 4, 0), C = c(3, 7, 1, 2, 1), E = c(5, 9, 1, 4, 0),
  Q = c(5, 10, 2, 3, 0), G = c(2, 5, 1, 2, 0), H = c(6, 9, 3, 2, 0),
  I = c(6, 13, 1, 2, 0), L = c(6, 13, 1, 2, 0), K = c(6, 14, 2, 2, 0),
  M = c(5, 11, 1, 2, 1), F = c(9, 11, 1, 2, 0), P = c(5, 9, 1, 2, 0),
  S = c(3, 7, 1, 3, 0), T = c(4, 9, 1, 3, 0), W = c(11, 12, 2, 2, 0),
  Y = c(9, 11, 1, 3, 0), V = c(5, 11, 1, 2, 0)
)

atomic_masses <- list(C = 12.000000, H = 1.007825, N = 14.003074, O = 15.994915, S = 31.972071)
proton_mass <- 1.007276

# 计算肽段分子式
calculate_peptide_formula <- function(sequence) {
  sequence <- toupper(sequence)
  aa_count <- table(strsplit(sequence, "")[[1]])
  total_formula <- c(C=0, H=0, N=0, O=0, S=0)
  
  for (aa in names(aa_count)) {
    if (aa %in% names(amino_acids)) {
      count <- aa_count[aa]
      formula <- amino_acids[[aa]]
      total_formula[1:5] <- total_formula[1:5] + formula * count
    }
  }
  
  water_loss <- nchar(sequence) - 1
  total_formula["H"] <- total_formula["H"] - water_loss * 2
  total_formula["O"] <- total_formula["O"] - water_loss * 1
  total_formula["H"] <- total_formula["H"] + 2
  total_formula["O"] <- total_formula["O"] + 1
  
  return(total_formula)
}

# 计算单同位素分子量
calculate_monoisotopic_mass <- function(formula) {
  mass <- 0
  for (element in names(formula)) {
    if (formula[element] > 0 && element %in% names(atomic_masses)) {
      mass <- mass + formula[element] * atomic_masses[[element]]
    }
  }
  return(mass)
}

# 二项分布概率计算
binomial_probability <- function(n, k, p) {
  if (k > n) return(0)
  if (k == 0) return((1-p)^n)
  if (k == n) return(p^n)
  
  log_result <- 0
  for (i in 1:k) {
    log_result <- log_result + log(n - i + 1) - log(i)
  }
  log_result <- log_result + k * log(p) + (n - k) * log(1 - p)
  return(exp(log_result))
}

# 计算同位素分布
calculate_isotope_distribution <- function(formula, charge = 1, top_n = 5) {
  mono_mass <- calculate_monoisotopic_mass(formula)
  nC <- formula["C"]
  pC13 <- 0.0107
  
  peaks <- data.frame()
  for (i in 0:min(10, nC)) {
    abundance <- binomial_probability(nC, i, pC13)
    isotope_mass <- mono_mass + i * 1.003355
    mz <- (isotope_mass + charge * proton_mass) / charge
    
    peaks <- rbind(peaks, data.frame(
      isotope = ifelse(i == 0, "M", paste0("M+", i)),
      neutral_mass = isotope_mass,
      mz = mz,
      abundance = abundance,
      relative_abundance = abundance * 100
    ))
  }
  
  peaks <- peaks[order(-peaks$abundance), ]
  return(peaks[1:min(top_n, nrow(peaks)), ])
}

# 计算综合分布 - 支持自然同位素和标记分析
calculate_comprehensive_distribution <- function(sequences, analysis_type = "natural", ratios = c(1), charge_range = 2:5) {
  all_results <- data.frame()
  
  for (seq in sequences) {
    seq <- trimws(seq)
    if (nchar(seq) == 0) next
    
    tryCatch({
      # 基础分子式
      formula <- calculate_peptide_formula(seq)
      base_mono_mass <- calculate_monoisotopic_mass(formula)
      
      # 创建分子式字符串
      formula_str <- paste0("C", formula["C"], "H", formula["H"], "N", formula["N"], 
                            "O", formula["O"], if(formula["S"] > 0) paste0("S", formula["S"]) else "")
      
      if (analysis_type == "natural") {
        # 自然同位素分析
        for (charge in charge_range) {
          # 计算自然同位素分布
          isotope_peaks <- calculate_isotope_distribution(formula, charge, 5)
          
          # 创建结果行
          row_data <- data.frame(
            Peptide_Sequence = seq,
            Length = nchar(seq),
            Charge_State = charge,
            Analysis_Type = "Natural",
            Molecular_Formula = formula_str,
            Monoisotopic_Mass = round(base_mono_mass, 4),
            Base_MZ = round((base_mono_mass + charge * proton_mass) / charge, 4),
            stringsAsFactors = FALSE
          )
          
          # 添加同位素峰数据 (M, M+1, M+2, M+3, M+4)
          for (i in 1:5) {
            isotope_col <- paste0("M", ifelse(i==1, "", paste0("_", i-1)))
            mz_col <- paste0("M", ifelse(i==1, "", paste0("_", i-1)), "_mz")
            
            if (i <= nrow(isotope_peaks)) {
              row_data[[isotope_col]] <- round(isotope_peaks$relative_abundance[i], 2)
              row_data[[mz_col]] <- round(isotope_peaks$mz[i], 4)
            } else {
              row_data[[isotope_col]] <- 0
              row_data[[mz_col]] <- 0
            }
          }
          
          all_results <- rbind(all_results, row_data)
        }
        
      } else {
        # 标记分析
        # 标记分子式 C18H25N3O3
        label_formula <- c(C=18, H=25, N=3, O=3, S=0)
        light_formula <- formula + label_formula
        heavy_formula <- light_formula  # 重标记只是质量差异
        
        for (charge in charge_range) {
          for (ratio in ratios) {
            # 计算轻重标记分布
            light_peaks <- calculate_isotope_distribution(light_formula, charge, 5)
            heavy_peaks <- calculate_isotope_distribution(heavy_formula, charge, 5)
            
            # 计算混合比例
            light_factor <- 1 / (1 + ratio)
            heavy_factor <- ratio / (1 + ratio)
            
            # 调整丰度
            light_peaks$adjusted_abundance <- light_peaks$abundance * light_factor * 100
            heavy_peaks$adjusted_abundance <- heavy_peaks$abundance * heavy_factor * 100
            
            # 计算质量信息
            light_mass <- base_mono_mass + 331.1896
            heavy_mass <- base_mono_mass + 333.1938
            light_base_mz <- (light_mass + charge * proton_mass) / charge
            heavy_base_mz <- (heavy_mass + charge * proton_mass) / charge
            
            # 创建结果行
            row_data <- data.frame(
              Peptide_Sequence = seq,
              Charge_State = charge,
              Analysis_Type = "Labeled",
              Label_Ratio = paste0("1:", ratio),
              Light_Percent = round(light_factor * 100, 1),
              Heavy_Percent = round(heavy_factor * 100, 1),
              Base_Mass = round(base_mono_mass, 4),
              Light_Mass = round(light_mass, 4),
              Heavy_Mass = round(heavy_mass, 4),
              Light_MZ = round(light_base_mz, 4),
              Heavy_MZ = round(heavy_base_mz, 4),
              stringsAsFactors = FALSE
            )
            
            # 添加轻标记峰 (M, M+1, M+2, M+3, M+4)
            for (i in 1:5) {
              if (i <= nrow(light_peaks)) {
                row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)))]] <- round(light_peaks$adjusted_abundance[i], 2)
                iso_mz <- (light_mass + (i-1)*1.003355 + charge * proton_mass) / charge
                row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- round(iso_mz, 4)
              } else {
                row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)))]] <- 0
                row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- 0
              }
            }
            
            # 添加重标记峰 (M', M'+1, M'+2, M'+3, M'+4)
            for (i in 1:5) {
              if (i <= nrow(heavy_peaks)) {
                row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)))]] <- round(heavy_peaks$adjusted_abundance[i], 2)
                iso_mz <- (heavy_mass + (i-1)*1.003355 + charge * proton_mass) / charge
                row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- round(iso_mz, 4)
              } else {
                row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)))]] <- 0
                row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- 0
              }
            }
            
            all_results <- rbind(all_results, row_data)
          }
        }
      }
    }, error = function(e) {
      warning(paste("Error processing sequence", seq, ":", e$message))
    })
    
    # 综合可视化 - 修复版本
    output$comprehensive_plot <- renderPlotly({
      req(comprehensive_data(), input$plot_peptide, input$plot_charge)
      
      data <- comprehensive_data()
      
      # 检查数据是否有Analysis_Type列
      if (!"Analysis_Type" %in% names(data)) {
        return(NULL)
      }
      
      analysis_type <- unique(data$Analysis_Type)[1]
      
      if (analysis_type == "Natural") {
        # 自然同位素可视化
        selected_data <- data[data$Peptide_Sequence == input$plot_peptide & 
                                data$Charge_State == as.numeric(input$plot_charge), ]
        
        if (nrow(selected_data) == 0) {
          return(plot_ly() %>% layout(title = "No data for selected conditions"))
        }
        
        # 准备可视化数据 - 使用实际的列名
        isotopes <- c("M", "M+1", "M+2", "M+3", "M+4")
        abundances <- c(
          if("M" %in% names(selected_data)) selected_data$M else 0,
          if("M_1" %in% names(selected_data)) selected_data$M_1 else 0,
          if("M_2" %in% names(selected_data)) selected_data$M_2 else 0,
          if("M_3" %in% names(selected_data)) selected_data$M_3 else 0,
          if("M_4" %in% names(selected_data)) selected_data$M_4 else 0
        )
        mz_values <- c(
          if("M_mz" %in% names(selected_data)) selected_data$M_mz else 0,
          if("M_1_mz" %in% names(selected_data)) selected_data$M_1_mz else 0,
          if("M_2_mz" %in% names(selected_data)) selected_data$M_2_mz else 0,
          if("M_3_mz" %in% names(selected_data)) selected_data$M_3_mz else 0,
          if("M_4_mz" %in% names(selected_data)) selected_data$M_4_mz else 0
        )
        
        # 创建数据框
        plot_data <- data.frame(
          Isotope = isotopes,
          Abundance = as.numeric(abundances),
          MZ = as.numeric(mz_values)
        )
        
        # 只保留丰度大于0.01%的峰
        plot_data <- plot_data[plot_data$Abundance > 0.01 & plot_data$MZ > 0, ]
        
        if (nrow(plot_data) == 0) {
          return(plot_ly() %>% layout(title = "No peaks found for visualization"))
        }
        
        p <- plot_ly(plot_data, x = ~MZ, y = ~Abundance, type = 'bar',
                     marker = list(color = 'steelblue', opacity = 0.8),
                     text = ~paste("Isotope:", Isotope, "<br>m/z:", MZ, "<br>Abundance:", round(Abundance, 2), "%"),
                     hovertemplate = "%{text}<extra></extra>") %>%
          layout(
            title = paste("Natural Isotope Distribution<br>",
                          input$plot_peptide, " | Charge: ", input$plot_charge, "+"),
            xaxis = list(title = "m/z"),
            yaxis = list(title = "Relative Abundance (%)"),
            showlegend = FALSE
          )
        
      } else {
        # 标记分析可视化
        req(input$plot_ratio)
        
        selected_data <- data[data$Peptide_Sequence == input$plot_peptide & 
                                data$Charge_State == as.numeric(input$plot_charge) & 
                                data$Label_Ratio == input$plot_ratio, ]
        
        if (nrow(selected_data) == 0) {
          return(plot_ly() %>% layout(title = "No data for selected conditions"))
        }
        
        # 准备轻重标记数据
        plot_data <- data.frame()
        
        # 轻标记数据
        light_isotopes <- c("M", "M+1", "M+2", "M+3", "M+4")
        light_abundances <- c(
          if("M" %in% names(selected_data)) selected_data$M else 0,
          if("M_1" %in% names(selected_data)) selected_data$M_1 else 0,
          if("M_2" %in% names(selected_data)) selected_data$M_2 else 0,
          if("M_3" %in% names(selected_data)) selected_data$M_3 else 0,
          if("M_4" %in% names(selected_data)) selected_data$M_4 else 0
        )
        light_mz <- c(
          if("M_mz" %in% names(selected_data)) selected_data$M_mz else 0,
          if("M_1_mz" %in% names(selected_data)) selected_data$M_1_mz else 0,
          if("M_2_mz" %in% names(selected_data)) selected_data$M_2_mz else 0,
          if("M_3_mz" %in% names(selected_data)) selected_data$M_3_mz else 0,
          if("M_4_mz" %in% names(selected_data)) selected_data$M_4_mz else 0
        )
        
        for (i in 1:5) {
          if (as.numeric(light_abundances[i]) > 0.01) {
            plot_data <- rbind(plot_data, data.frame(
              Isotope = light_isotopes[i],
              Abundance = as.numeric(light_abundances[i]),
              MZ = as.numeric(light_mz[i]),
              Label = "Light"
            ))
          }
        }
        
        # 重标记数据
        heavy_isotopes <- c("M'", "M'+1", "M'+2", "M'+3", "M'+4")
        heavy_abundances <- c(
          if("Mp" %in% names(selected_data)) selected_data$Mp else 0,
          if("Mp_1" %in% names(selected_data)) selected_data$Mp_1 else 0,
          if("Mp_2" %in% names(selected_data)) selected_data$Mp_2 else 0,
          if("Mp_3" %in% names(selected_data)) selected_data$Mp_3 else 0,
          if("Mp_4" %in% names(selected_data)) selected_data$Mp_4 else 0
        )
        heavy_mz <- c(
          if("Mp_mz" %in% names(selected_data)) selected_data$Mp_mz else 0,
          if("Mp_1_mz" %in% names(selected_data)) selected_data$Mp_1_mz else 0,
          if("Mp_2_mz" %in% names(selected_data)) selected_data$Mp_2_mz else 0,
          if("Mp_3_mz" %in% names(selected_data)) selected_data$Mp_3_mz else 0,
          if("Mp_4_mz" %in% names(selected_data)) selected_data$Mp_4_mz else 0
        )
        
        for (i in 1:5) {
          if (as.numeric(heavy_abundances[i]) > 0.01) {
            plot_data <- rbind(plot_data, data.frame(
              Isotope = heavy_isotopes[i],
              Abundance = as.numeric(heavy_abundances[i]),
              MZ = as.numeric(heavy_mz[i]),
              Label = "Heavy"
            ))
          }
        }
        
        if (nrow(plot_data) == 0) {
          return(plot_ly() %>% layout(title = "No peaks found for visualization"))
        }
        
        p <- plot_ly(plot_data, x = ~MZ, y = ~Abundance, 
                     color = ~Label, colors = c('lightblue', 'salmon'),
                     type = 'bar',
                     text = ~paste("Isotope:", Isotope, "<br>m/z:", MZ, "<br>Abundance:", round(Abundance, 2), "%"),
                     hovertemplate = "%{text}<extra></extra>") %>%
          layout(
            title = paste("Labeled Peptide Distribution<br>",
                          input$plot_peptide, " | Charge: ", input$plot_charge, "+ | Ratio: ", input$plot_ratio),
            xaxis = list(title = "m/z"),
            yaxis = list(title = "Relative Abundance (%)"),
            barmode = 'group'
          )
      }
      
      return(p)
    })
    
    # 自然同位素CSV下载
    output$download_natural_csv <- downloadHandler(
      filename = function() {
        paste("natural_isotope_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
      },
      content = function(file) {
        req(natural_data())
        write.csv(natural_data(), file, row.names = FALSE)
      }
    )
  }
  
  return(all_results)
}
calculate_natural_distribution <- function(sequences, charge_range = 2:5) {
  all_results <- data.frame()
  
  for (seq in sequences) {
    seq <- trimws(seq)
    if (nchar(seq) == 0) next
    
    tryCatch({
      # 计算肽段分子式
      formula <- calculate_peptide_formula(seq)
      mono_mass <- calculate_monoisotopic_mass(formula)
      
      # 创建分子式字符串
      formula_str <- paste0("C", formula["C"], "H", formula["H"], "N", formula["N"], 
                            "O", formula["O"], if(formula["S"] > 0) paste0("S", formula["S"]) else "")
      
      for (charge in charge_range) {
        # 计算同位素分布
        isotope_peaks <- calculate_isotope_distribution(formula, charge, 5)
        
        # 创建结果行
        row_data <- data.frame(
          Peptide_Sequence = seq,
          Length = nchar(seq),
          Charge_State = charge,
          Molecular_Formula = formula_str,
          Monoisotopic_Mass = round(mono_mass, 4),
          Base_MZ = round((mono_mass + charge * proton_mass) / charge, 4),
          stringsAsFactors = FALSE
        )
        
        # 添加同位素峰数据 (M, M+1, M+2, M+3, M+4)
        for (i in 1:5) {
          isotope_col <- paste0("M", ifelse(i==1, "", paste0("_", i-1)))
          mz_col <- paste0("M", ifelse(i==1, "", paste0("_", i-1)), "_mz")
          
          if (i <= nrow(isotope_peaks)) {
            row_data[[isotope_col]] <- round(isotope_peaks$relative_abundance[i], 2)
            row_data[[mz_col]] <- round(isotope_peaks$mz[i], 4)
          } else {
            row_data[[isotope_col]] <- 0
            row_data[[mz_col]] <- 0
          }
        }
        
        all_results <- rbind(all_results, row_data)
      }
    }, error = function(e) {
      warning(paste("Error processing sequence", seq, ":", e$message))
    })
    
    # 自然同位素CSV下载
    output$download_natural_csv <- downloadHandler(
      filename = function() {
        paste("natural_isotope_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
      },
      content = function(file) {
        req(natural_data())
        write.csv(natural_data(), file, row.names = FALSE)
      }
    )
    
    # 自然同位素表格
    output$natural_table <- DT::renderDataTable({
      req(natural_data())
      data <- natural_data()
      
      DT::datatable(data, 
                    options = list(pageLength = 15, scrollX = TRUE, scrollY = "400px"),
                    rownames = FALSE) %>%
        # 同位素丰度列 - 蓝色背景
        DT::formatStyle(columns = c("M", "M_1", "M_2", "M_3", "M_4"),
                        backgroundColor = 'rgba(135, 206, 250, 0.3)',
                        fontWeight = 'bold') %>%
        # m/z列 - 绿色背景
        DT::formatStyle(columns = grep("_mz$", names(data), value = TRUE),
                        backgroundColor = 'rgba(152, 251, 152, 0.3)',
                        fontSize = '12px') %>%
        # 基础信息列 - 浅灰色背景
        DT::formatStyle(columns = c("Molecular_Formula", "Monoisotopic_Mass", "Base_MZ"),
                        backgroundColor = 'rgba(245, 245, 245, 0.8)')
    })
    
    # 自然同位素摘要
    output$natural_summary <- renderText({
      if (is.null(natural_data())) {
        return("No analysis results yet. Enter sequences and click 'Calculate'.")
      }
      
      data <- natural_data()
      unique_sequences <- unique(data$Peptide_Sequence)
      total_combinations <- nrow(data)
      
      paste(
        "Analysis Summary:",
        paste("• Peptides processed:", length(unique_sequences)),
        paste("• Charge states:", paste(sort(unique(data$Charge_State)), collapse = ", ")),
        paste("• Total combinations:", total_combinations),
        paste("• Natural isotope peaks: M to M+4"),
        sep = "\n"
      )
    })
    
    # 自然同位素可视化
    output$natural_plot <- renderPlotly({
      req(natural_data(), input$plot_natural_peptide, input$plot_natural_charge)
      
      data <- natural_data()
      
      # 筛选选中条件的数据
      selected_data <- data[data$Peptide_Sequence == input$plot_natural_peptide & 
                              data$Charge_State == as.numeric(input$plot_natural_charge), ]
      
      if (nrow(selected_data) == 0) return(NULL)
      
      # 准备可视化数据
      plot_data <- data.frame()
      
      isotopes <- c("M", "M+1", "M+2", "M+3", "M+4")
      abundances <- c(selected_data$M, selected_data$M_1, selected_data$M_2, 
                      selected_data$M_3, selected_data$M_4)
      mz_values <- c(selected_data$M_mz, selected_data$M_1_mz, selected_data$M_2_mz,
                     selected_data$M_3_mz, selected_data$M_4_mz)
      
      for (i in 1:5) {
        if (abundances[i] > 0.01) {  # 只显示丰度大于0.01%的峰
          plot_data <- rbind(plot_data, data.frame(
            Isotope = isotopes[i],
            Abundance = abundances[i],
            MZ = mz_values[i]
          ))
        }
      }
      
      # 创建图表
      p <- plot_ly(plot_data, x = ~MZ, y = ~Abundance, type = 'bar',
                   marker = list(color = 'steelblue', opacity = 0.8),
                   text = ~paste("Isotope:", Isotope, "<br>m/z:", MZ, "<br>Abundance:", round(Abundance, 2), "%"),
                   hovertemplate = "%{text}<extra></extra>") %>%
        layout(
          title = paste("Natural Isotope Distribution", "<br>",
                        input$plot_natural_peptide, "| Charge:", paste0(input$plot_natural_charge, "+"), 
                        "| Formula:", selected_data$Molecular_Formula),
          xaxis = list(title = "m/z", showgrid = TRUE),
          yaxis = list(title = "Relative Abundance (%)", showgrid = TRUE),
          plot_bgcolor = 'rgba(240,240,240,0.3)',
          showlegend = FALSE
        )
      
      return(p)
    })
    
    # 自然同位素分析
    observeEvent(input$calculate_natural, {
      req(input$natural_charges, input$natural_sequences)
      
      sequences <- strsplit(input$natural_sequences, "\n")[[1]]
      sequences <- sequences[nchar(trimws(sequences)) > 0]
      
      if (length(sequences) == 0) {
        showNotification("Please provide sequences", type = "warning")
        return()
      }
      
      showNotification("Calculating natural isotopes...", type = "message")
      
      results <- calculate_natural_distribution(sequences, as.numeric(input$natural_charges))
      natural_data(results)
      
      # 更新可视化选择器
      unique_peptides <- unique(results$Peptide_Sequence)
      unique_charges <- sort(unique(results$Charge_State))
      
      updateSelectInput(session, "plot_natural_peptide", 
                        choices = setNames(unique_peptides, unique_peptides),
                        selected = unique_peptides[1])
      updateSelectInput(session, "plot_natural_charge", 
                        choices = setNames(unique_charges, paste0(unique_charges, "+")),
                        selected = unique_charges[1])
      
      showNotification("Natural isotope analysis completed!", type = "message")
    })
  }
  
  return(all_results)
}
calculate_comprehensive_distribution <- function(sequences, ratios = c(1), charge_range = 2:5) {
  all_results <- data.frame()
  
  for (seq in sequences) {
    seq <- trimws(seq)
    if (nchar(seq) == 0) next
    
    tryCatch({
      # 基础分子式
      formula <- calculate_peptide_formula(seq)
      base_mono_mass <- calculate_monoisotopic_mass(formula)
      
      # 标记分子式 C18H25N3O3
      label_formula <- c(C=18, H=25, N=3, O=3, S=0)
      light_formula <- formula + label_formula
      heavy_formula <- light_formula  # 重标记只是质量差异
      
      for (charge in charge_range) {
        for (ratio in ratios) {
          # 计算轻重标记分布
          light_peaks <- calculate_isotope_distribution(light_formula, charge, 5)
          heavy_peaks <- calculate_isotope_distribution(heavy_formula, charge, 5)
          
          # 计算混合比例
          light_factor <- 1 / (1 + ratio)
          heavy_factor <- ratio / (1 + ratio)
          
          # 调整丰度
          light_peaks$adjusted_abundance <- light_peaks$abundance * light_factor * 100
          heavy_peaks$adjusted_abundance <- heavy_peaks$abundance * heavy_factor * 100
          
          # 计算质量信息
          light_mass <- base_mono_mass + 331.1896
          heavy_mass <- base_mono_mass + 333.1938
          light_base_mz <- (light_mass + charge * proton_mass) / charge
          heavy_base_mz <- (heavy_mass + charge * proton_mass) / charge
          
          # 创建结果行
          row_data <- data.frame(
            Peptide_Sequence = seq,
            Charge_State = charge,
            Label_Ratio = paste0("1:", ratio),
            Light_Percent = round(light_factor * 100, 1),
            Heavy_Percent = round(heavy_factor * 100, 1),
            Base_Mass = round(base_mono_mass, 4),
            Light_Mass = round(light_mass, 4),
            Heavy_Mass = round(heavy_mass, 4),
            Light_MZ = round(light_base_mz, 4),
            Heavy_MZ = round(heavy_base_mz, 4),
            stringsAsFactors = FALSE
          )
          
          # 添加轻标记峰 (M, M+1, M+2, M+3, M+4)
          for (i in 1:5) {
            if (i <= nrow(light_peaks)) {
              row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)))]] <- round(light_peaks$adjusted_abundance[i], 2)
              iso_mz <- (light_mass + (i-1)*1.003355 + charge * proton_mass) / charge
              row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- round(iso_mz, 4)
            } else {
              row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)))]] <- 0
              row_data[[paste0("M", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- 0
            }
          }
          
          # 添加重标记峰 (M', M'+1, M'+2, M'+3, M'+4)
          for (i in 1:5) {
            if (i <= nrow(heavy_peaks)) {
              row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)))]] <- round(heavy_peaks$adjusted_abundance[i], 2)
              iso_mz <- (heavy_mass + (i-1)*1.003355 + charge * proton_mass) / charge
              row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- round(iso_mz, 4)
            } else {
              row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)))]] <- 0
              row_data[[paste0("Mp", ifelse(i==1, "", paste0("_", i-1)), "_mz")]] <- 0
            }
          }
          
          all_results <- rbind(all_results, row_data)
        }
      }
    }, error = function(e) {
      warning(paste("Error processing sequence", seq, ":", e$message))
    })
  }
  
  return(all_results)
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Peptide Analyzer"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Basic Analysis", tabName = "basic", icon = icon("calculator")),
      menuItem("Comprehensive Analysis", tabName = "comprehensive", icon = icon("table"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # 基础分析
      tabItem(tabName = "basic",
              fluidRow(
                box(
                  title = "Input", status = "primary", solidHeader = TRUE, width = 4,
                  textInput("sequence", "Peptide Sequence:", value = "HIADLAGNSEVILPVPAFNVINGGSHAGNK"),
                  numericInput("charge", "Charge State:", value = 3, min = 1, max = 10),
                  actionButton("calculate", "Calculate", class = "btn-primary")
                ),
                
                # 自然同位素分析
                tabItem(tabName = "natural",
                        fluidRow(
                          box(
                            title = "Natural Isotope Analysis", status = "primary", solidHeader = TRUE, width = 4,
                            p("Calculate natural isotope distribution without any labels"),
                            textAreaInput("natural_sequences", "Peptide Sequences (one per line):", 
                                          value = "HIADLAGNSEVILPVPAFNVINGGSHAGNK\nASPILGIDVWEIHELLAK\nPEPTIDE\nPROTEIN", 
                                          rows = 8),
                            checkboxGroupInput("natural_charges", "Charge States:", 
                                               choices = list("2+" = 2, "3+" = 3, "4+" = 4, "5+" = 5),
                                               selected = c(2, 3, 4, 5), inline = TRUE),
                            actionButton("calculate_natural", "Calculate Natural Isotopes", class = "btn-success"),
                            hr(),
                            verbatimTextOutput("natural_summary"),
                            downloadButton("download_natural_csv", "Download CSV", class = "btn-info")
                          ),
                          box(
                            title = "Natural Isotope Results", status = "info", solidHeader = TRUE, width = 8,
                            p("Natural isotope distribution showing M, M+1, M+2, M+3, M+4 peaks with corresponding m/z values"),
                            tabsetPanel(
                              tabPanel("Data Table",
                                       DT::dataTableOutput("natural_table")
                              ),
                              tabPanel("Visualization",
                                       fluidRow(
                                         column(6, selectInput("plot_natural_peptide", "Select Peptide:", choices = NULL)),
                                         column(6, selectInput("plot_natural_charge", "Select Charge:", choices = NULL))
                                       ),
                                       plotlyOutput("natural_plot", height = "500px")
                              )
                            )
                          )
                        )
                ),
                box(
                  title = "Results", status = "success", solidHeader = TRUE, width = 8,
                  DT::dataTableOutput("basic_table"),
                  plotlyOutput("basic_plot")
                )
              )
      ),
      
      # 综合分析
      tabItem(tabName = "comprehensive",
              fluidRow(
                box(
                  title = "Comprehensive Analysis", status = "primary", solidHeader = TRUE, width = 4,
                  p("Select analysis type:"),
                  radioButtons("analysis_type", "Analysis Type:",
                               choices = list(
                                 "Natural Isotopes (no labels)" = "natural",
                                 "Labeled Peptides (light/heavy)" = "labeled"
                               ),
                               selected = "natural"),
                  hr(),
                  
                  textAreaInput("batch_sequences", "Peptide Sequences (one per line):", 
                                value = "HIADLAGNSEVILPVPAFNVINGGSHAGNK\nASPILGIDVWEIHELLAK\nPEPTIDE\nPROTEIN", 
                                rows = 6),
                  
                  checkboxGroupInput("charge_states", "Charge States:", 
                                     choices = list("2+" = 2, "3+" = 3, "4+" = 4, "5+" = 5),
                                     selected = c(2, 3), inline = TRUE),
                  
                  # 条件显示：只有选择labeled时才显示比例选项
                  conditionalPanel(
                    condition = "input.analysis_type == 'labeled'",
                    hr(),
                    p(strong("Label Ratios (Light:Heavy):")),
                    checkboxGroupInput("ratios", "Ratios:", 
                                       choices = list("1:1" = 1, "1:2" = 2, "1:3" = 3, "1:4" = 4, "1:5" = 5, "1:6" = 6),
                                       selected = c(1, 6), inline = TRUE)
                  ),
                  
                  hr(),
                  actionButton("calculate_comprehensive", "Calculate Analysis", class = "btn-primary"            ),
                  hr(),
                  verbatimTextOutput("analysis_summary"),
                  br(),
                  div(style = "font-size: 12px; color: gray;",
                      textOutput("debug_info")
                  ),
                  downloadButton("download_csv", "Download CSV", class = "btn-success")
                ),
                
                box(
                  title = "Analysis Results", status = "warning", solidHeader = TRUE, width = 8,
                  
                  # 动态显示说明文字
                  uiOutput("results_description"),
                  
                  tabsetPanel(
                    tabPanel("Data Table",
                             DT::dataTableOutput("comprehensive_table")
                    ),
                    tabPanel("Visualization",
                             # 动态UI根据分析类型显示不同的选择器
                             uiOutput("plot_controls"),
                             plotlyOutput("comprehensive_plot", height = "500px")
                    )
                  )
                )
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  basic_data <- reactiveVal()
  comprehensive_data <- reactiveVal()
  
  # 基础分析
  observeEvent(input$calculate, {
    req(input$sequence, input$charge)
    formula <- calculate_peptide_formula(input$sequence)
    peaks <- calculate_isotope_distribution(formula, input$charge, 5)
    basic_data(peaks)
  })
  
  # 综合分析
  observeEvent(input$calculate_comprehensive, {
    req(input$charge_states, input$ratios, input$batch_sequences)
    
    sequences <- strsplit(input$batch_sequences, "\n")[[1]]
    sequences <- sequences[nchar(trimws(sequences)) > 0]
    
    if (length(sequences) == 0) {
      showNotification("Please provide sequences", type = "warning")
      return()
    }
    
    showNotification("Processing comprehensive analysis...", type = "message")
    
    results <- calculate_comprehensive_distribution(sequences, as.numeric(input$ratios), as.numeric(input$charge_states))
    comprehensive_data(results)
    
    showNotification("Comprehensive analysis completed!", type = "message")
  })
  
  # 基础表格
  output$basic_table <- DT::renderDataTable({
    req(basic_data())
    data <- basic_data()
    data$neutral_mass <- round(data$neutral_mass, 4)
    data$mz <- round(data$mz, 4)
    data$relative_abundance <- round(data$relative_abundance, 2)
    DT::datatable(data[, c("isotope", "mz", "relative_abundance")], 
                  options = list(dom = 't'), rownames = FALSE)
  })
  
  # 综合表格
  output$comprehensive_table <- DT::renderDataTable({
    req(comprehensive_data())
    data <- comprehensive_data()
    
    DT::datatable(data, 
                  options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE) %>%
      DT::formatStyle(columns = c("M", "M_1", "M_2", "M_3", "M_4"),
                      backgroundColor = 'rgba(173, 216, 230, 0.3)') %>%
      DT::formatStyle(columns = c("Mp", "Mp_1", "Mp_2", "Mp_3", "Mp_4"),
                      backgroundColor = 'rgba(250, 128, 114, 0.3)') %>%
      DT::formatStyle(columns = grep("_mz$", names(data), value = TRUE),
                      backgroundColor = 'rgba(144, 238, 144, 0.2)')
  })
  
  # 基础图表
  output$basic_plot <- renderPlotly({
    req(basic_data())
    data <- basic_data()
    
    plot_ly(data, x = ~mz, y = ~relative_abundance, type = 'bar',
            text = ~isotope, textposition = 'outside') %>%
      layout(title = "Isotope Distribution", 
             xaxis = list(title = "m/z"),
             yaxis = list(title = "Abundance (%)"))
  })
  
  # CSV下载
  output$download_csv <- downloadHandler(
    filename = function() {
      paste("comprehensive_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
    },
    content = function(file) {
      req(comprehensive_data())
      write.csv(comprehensive_data(), file, row.names = FALSE)
    }
  )
}

# 运行应用
shinyApp(ui = ui, server = server)

shinyApp(ui = ui, server = server)