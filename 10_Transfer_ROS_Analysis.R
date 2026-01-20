# ==============================================================================
# Title: Transfer & ROS Analysis (P1: Pesticide diversity groups | P2: Pesticide diversity levels)
# Description: Profiling with adaptive jitter mapping, 
#              enlarged data points, and statistical significance testing.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

# --- Define Professional Color Palettes ---

# D0-D4 Colors (Specified RGB)
group2_colors <- c(
  "D0" = "#B11AB5",
  "D1" = "#EECA40",
  "D2" = "#FD763F",
  "D3" = "#23BAC5",
  "D4" = "#65BA45"
  )

# CK, Low, High Diversity Treatment Colors (For P1 Bars)
treatment_palette <- c(
  "CK"   = "#7B9CEE", 
  "Low"  = "#E88571", 
  "High" = "#82B461"  
)

# Mixture Subgroup Colors (For P2 Jitter)
mixture_palette <- c(
  "A"="#C68BB0", "AB"="#E9C487", "ABC"="#96C8CD", "ABCD"="#A3CC8B", "ABD"="#95C6BD",
  "AC"="#E5B17E", "ACD"="#97C7A8", "AD"="#E1A17B", "B"="#CE9E9C", "BC"="#D4A184",
  "BCD"="#9CC898", "BD"="#B9B09D", "C"="#DBB98D", "CD"="#A2B9B3", "CK"="#B16BBF", "D"="#E9D28B"
)

# 2. Data Loading and Matrix Preprocessing
# ------------------------------------------------------------------------------
# - Columns 1-4: Metadata (Column 1: sample ID, Column 2: group1 (CK,Low,High), Column 3: group2 (D0-4), Column 4: Mixture (specific pesticide combination)
# - Column 5+: Features (Numeric values for Transfer & ROS)
# - Rows: Each row represents a unique Sample coupled with its Grouping info.
df <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

target_columns <- c("Transfer frequency (x104)", 
                    "Fluorescence intensity on ROS levels of donor", 
                    "Fluorescence intensity on ROS levels of recipient")

# Global Factor Level Enforcement
if("group1" %in% colnames(df)) df$group1 <- factor(df$group1, levels = c("CK", "Low", "High"))
if("group" %in% colnames(df))  df$group2 <- factor(df$group2,  levels = c("D0", "D1", "D2", "D3", "D4"))

# 3. Core Visualization Function with Adaptive Mapping
# ------------------------------------------------------------------------------
draw_adaptive_profile_plot <- function(x_axis_selection = "group2") {
  
  # A. Conditional logic for jitter and fill mapping
  color_var    <- ifelse(x_axis_selection == "group1", "group2", "Mixture")
  color_scales <- if(color_var == "group2") group2_colors else mixture_palette
  fill_scales  <- if(x_axis_selection == "group1") treatment_palette else group2_colors
  
  # B. Processing Data for Multi-Panel Faceting
  df_long <- df %>%
    pivot_longer(cols = all_of(target_columns), names_to = "Metric", values_to = "Value") %>%
    mutate(metric_label = str_wrap(Metric, width = 25))
  
  # C. Statistical Calculations for Bars and Error Bars
  df_bar_stats <- df_long %>%
    group_by(!!sym(x_axis_selection), metric_label) %>%
    summarise(mean_val = mean(Value, na.rm = TRUE),
              se_val   = sd(Value, na.rm = TRUE) / sqrt(n()), .groups = "drop")
  
  # D. Plotting
  p <- ggplot(df_long, aes(x = !!sym(x_axis_selection), y = Value)) +
    facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
    
    # Standard Bars (Mean)
    geom_col(data = df_bar_stats, aes(y = mean_val, fill = !!sym(x_axis_selection)), 
             width = 0.65, color = "black", size = 0.5, alpha = 0.8) +
    
    # Error Bars (SE)
    geom_errorbar(data = df_bar_stats, 
                  aes(y = mean_val, ymin = mean_val - se_val, ymax = mean_val + se_val), 
                  width = 0.18, size = 0.7, color = "black") +
    
    # ENLARGED Jitter Points
    geom_jitter(aes(color = !!sym(color_var)), 
                width = 0.2, size = 3.8, alpha = 0.8, stroke = 0.4) +
    
    # 
    
    scale_fill_manual(values = fill_scales) +
    scale_color_manual(values = color_scales) +
    theme_classic() +
    theme(
      axis.text   = element_text(size = 11, color = "black"),
      axis.title  = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "gray92", color = "black"),
      strip.text  = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.spacing = unit(1.2, "lines")
    ) +
    labs(
      x = NULL, 
      y = "Measured Values (Mean ± SE)",
      fill = ifelse(x_axis_selection == "group1", "Treatment", "Mixture"),
      color = color_var,
      title = paste("Profile Analysis by", x_axis_selection)
    )
  
  # E. Statistics: Adaptive Comparisons
  if(x_axis_selection == "group1") {
    # P1: Comparisons for Diversity groups
    p <- p + stat_compare_means(method = "t.test", 
                                comparisons = list(c("CK", "Low"), c("CK", "High")),
                                label = "p.signif")
  } else {
    # P2: Comparisons for Diversity levels (D0-D2)
    p <- p + stat_compare_means(data = df_long %>% filter(group2 %in% c("D0", "D1", "D2")),
                                method = "wilcox.test", 
                                comparisons = list(c("D0", "D1"), c("D1", "D2"), c("D0", "D2")),
                                label = "p.signif", step.increase = 0.12)
  }
  
  return(p)
}

# 4. Final Output Generation
# ------------------------------------------------------------------------------
# P1: Focused on Pesticide Diversity groups (CK, Low, High)
p_pesticide_groups <- draw_adaptive_profile_plot(x_axis_selection = "group1")

# P2: Focused on Pesticide Diversity levels (D0, D1, D2, D3, D4)
p_pesticide_levels  <- draw_adaptive_profile_plot(x_axis_selection = "group2")

# Display Results
print(p_pesticide_groups)
print(p_pesticide_levels)
