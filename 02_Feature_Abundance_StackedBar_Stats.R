# ==============================================================================
# Title: Feature Profile Analysis (Stacked Bar + Jitter + T-test)
# Description: Profiling of features (ARG/MGE/VFG) with adaptive 
#              jitter mapping and statistical significance testing.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

# Define color palettes
feature_palette <- c("#7F95D1", "#57D1C9", "#FDC4B6", "#D9534F", "#CAB2D6", 
                     "#2694AB", "#B3DE69", "#96CEB4", "#FFEEAD", "#EA7070", "#8DD3C7")

group2_colors <- c(
  "D0" = "#7B81AF", "D1" = "#D3D3D3", "D2" = "#73AB81", 
  "D3" = "#87619D", "D4" = "#72A6A5"
)

#color palette for active metagenome
#group2_colors <- c(
#  "D0" = "#B11AB5",
#  "D1" = "#EECA40",
#  "D2" = "#FD763F",
#  "D3" = "#23BAC5",
#  "D4" = "#65BA45"
#)

# 2. Data Loading and Matrix Preprocessing
# ------------------------------------------------------------------------------
# - Columns 1-4: Metadata (Column 1: sample ID, Column 2: group1 (CK,Low,High), Column 3: group2 (D0-4), Column 4: Mixture (specific pesticide combination)
# - Column 5+: Features (Numeric values for Taxa, Genes, or MGE types)
# - Rows: Each row represents a unique Sample coupled with its Grouping info.
df <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

val_cols <- df %>% select(where(is.numeric)) %>% colnames()

if("group1" %in% colnames(df)) df$group1 <- factor(df$group1, levels = c("CK", "Low", "High"))
if("group2" %in% colnames(df)) df$group2 <- factor(df$group2, levels = c("D0", "D1", "D2", "D3", "D4"))

# 3. Core Visualization Function with Adaptive Mapping
# ------------------------------------------------------------------------------
draw_adaptive_stacked_plot <- function(x_axis_selection = "group1", feature_label = "ARG") {
  
  # A. Conditional logic for jitter coloring
  color_var <- ifelse(x_axis_selection == "group1", "group2", "Mixture")
  
  # B. Processing Stacked Bar Data (Top 9 Filtering)
  df_long <- df %>%
    pivot_longer(cols = all_of(val_cols), names_to = "FeatureClass", values_to = "Abundance")
  
  top9 <- df_long %>%
    group_by(FeatureClass) %>% 
    summarise(m = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>% slice_head(n = 9) %>% pull(FeatureClass)
  
  df_bar_stats <- df_long %>%
    mutate(FeatureClass = ifelse(FeatureClass %in% top9, FeatureClass, "Others")) %>%
    group_by(!!sym(x_axis_selection), FeatureClass) %>%
    summarise(mean_val = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  # C. Statistical Calculations for Jitter and Error Bars
  df_sample_totals <- df %>%
    mutate(total_sum = rowSums(select(., all_of(val_cols)), na.rm = TRUE))
  
  df_error_stats <- df_bar_stats %>%
    group_by(!!sym(x_axis_selection)) %>%
    summarise(total_height = sum(mean_val), .groups = "drop") %>%
    left_join(
      df_sample_totals %>% group_by(!!sym(x_axis_selection)) %>%
        summarise(se = sd(total_sum, na.rm = TRUE) / sqrt(n()), .groups = "drop"),
      by = x_axis_selection
    )
  
  # D. Plotting
  y_limit <- max(df_sample_totals$total_sum, df_error_stats$total_height + df_error_stats$se) * 1.15
  
  p <- ggplot(df_sample_totals, aes(x = !!sym(x_axis_selection), y = total_sum)) +
    # Stacked Bars
    geom_col(data = df_bar_stats, aes(y = mean_val, fill = FeatureClass), width = 0.6) +
    # Error Bars - FIXED: Added explicit x mapping
    geom_errorbar(data = df_error_stats, 
                  aes(x = !!sym(x_axis_selection), 
                      ymin = total_height - se, 
                      ymax = total_height + se), 
                  width = 0.12, size = 0.7, color = "black", inherit.aes = FALSE) +
    # Jitter Points
    geom_jitter(aes(color = !!sym(color_var)), width = 0.15, size = 2.5, alpha = 0.8) +
    scale_fill_manual(values = feature_palette) +
    theme_classic() +
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(
      x = NULL, 
      y = paste(feature_label, "Abundance (Mean ± SE)"),
      fill = paste(feature_label, "Class"),
      color = color_var,
      title = paste("Abundance Profile by", x_axis_selection)
    )
  
  # E. Statistics: Add only if x_axis is group1
  if(x_axis_selection == "group1") {
    p <- p + stat_compare_means(
      method = "t.test",
      comparisons = list(c("CK", "Low"), c("CK", "High")),
      label = "p.signif",
      label.y = y_limit
    )
  }
  
  # Apply professional group2 colors if mapping to group2
  if(color_var == "group2") {
    p <- p + scale_color_manual(values = group2_colors)
  }
  
  return(p)
}

# 4. Final Output
# ------------------------------------------------------------------------------
p_group1 <- draw_adaptive_stacked_plot(x_axis_selection = "group1")
p_group2 <- draw_adaptive_stacked_plot(x_axis_selection = "group2")

print(p_group1)
print(p_group2)