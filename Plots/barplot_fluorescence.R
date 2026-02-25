#!/usr/bin/env Rscript
# ==============================================================================
# Bar Plot - Fluorescence Intensity (Duplicates with Error Bars)
# Phenol vs Azide distinguished by bar border linetype (solid vs dashed)
# FITC, EY, FCPP, ECPP distinguished by fill color
# ==============================================================================

# Load libraries
if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cran.r-project.org")
if (!require("dplyr")) install.packages("dplyr", repos = "https://cran.r-project.org")

library(ggplot2)
library(dplyr)

# ---- Raw data (duplicates) ----
raw_data <- data.frame(
  Sample = c(
    "FITC_Phenol", "FITC_Phenol",
    "EY_Phenol",   "EY_Phenol",
    "FITC_Azide",  "FITC_Azide",
    "EY_Azide",    "EY_Azide",
    "FCPP_Phenol", "FCPP_Phenol",
    "ECPP_Phenol", "ECPP_Phenol",
    "FCPP_Azide",  "FCPP_Azide",
    "ECPP_Azide",  "ECPP_Azide"
  ),
  Value = c(
    975, 1055,
    1554, 1532,
    995, 1034,
    2074, 2551,
    1646, 1662,
    247, 90,
    658, 626,
    632, 432
  )
)

# ---- Parse sample names into Dye and Condition ----
raw_data <- raw_data %>%
  mutate(
    Dye       = sub("_.*", "", Sample),
    Condition = sub(".*_", "", Sample)
  )

# ---- Summarize: mean and SD for each group ----
plot_data <- raw_data %>%
  group_by(Dye, Condition) %>%
  summarise(
    Mean = mean(Value),
    SD   = sd(Value),
    .groups = "drop"
  )

# ---- Set factor levels for desired plotting order ----
plot_data$Dye <- factor(plot_data$Dye, levels = c("FITC", "EY", "FCPP", "ECPP"))
plot_data$Condition <- factor(plot_data$Condition, levels = c("Phenol", "Azide"))

# ---- Define colors for each dye ----
dye_colors <- c(
  "FITC" = "#4DAF4A",   # green
  "EY"   = "#FF7F00",   # orange
  "FCPP" = "#377EB8",   # blue
  "ECPP" = "#E41A1C"    # red
)

# ---- Build the plot ----
p <- ggplot(plot_data, aes(x = Dye, y = Mean, fill = Dye, linetype = Condition)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.7,
    color    = "black",   # bar border color
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.8),
    width    = 0.25,
    linewidth = 0.6,
    color    = "black"
  ) +
  scale_fill_manual(values = dye_colors) +
  scale_linetype_manual(
    values = c("Phenol" = "solid", "Azide" = "dashed"),
    name   = "Condition"
  ) +
  labs(
    x     = "",
    y     = "Fluorescence Intensity (a.u.)",
    fill  = "Dye",
    title = "Fluorescence Intensity by Dye and Condition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    legend.position   = "right",
    axis.text.x       = element_text(size = 12, face = "bold"),
    axis.text.y       = element_text(size = 11),
    axis.title.y      = element_text(size = 13)
  ) +
  # Override the fill legend to also show linetype in the guide
  guides(
    fill     = guide_legend(override.aes = list(linetype = 0)),
    linetype = guide_legend(override.aes = list(fill = "grey80"))
  )

# ---- Save to PDF and PNG ----
output_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# Fallback: if not running in RStudio, save to script directory
if (is.null(output_dir) || output_dir == "") {
  output_dir <- "."
}

ggsave(
  filename = file.path(".", "barplot_fluorescence.pdf"),
  plot     = p,
  width    = 8,
  height   = 5,
  dpi      = 300
)

ggsave(
  filename = file.path(".", "barplot_fluorescence.png"),
  plot     = p,
  width    = 8,
  height   = 5,
  dpi      = 300
)

cat("Plot saved as barplot_fluorescence.pdf and barplot_fluorescence.png\n")
