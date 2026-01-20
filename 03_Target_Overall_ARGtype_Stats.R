# ==============================================================================
# Title: Targeted overall ARGtype Comparison (Violin + Boxplot + Outlier Labeling)
# Description: Detailed abundance analysis of specific genes using t-tests, 
#              violin plots, and automated outlier detection/labeling.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(ggrepel)

# Define professional RGB color palette
group2_colors <- c(
  "D0"      = "#7B81AF", # Baseline
  "D1"      = "#D3D3D3", 
  "D2"      = "#73AB81", 
  "D3"      = "#87619D", 
  "D4"      = "#72A6A5", # End point
  "Low"     = "#E0E0E0", # Consolidated Group
  "Outlier" = "#808080"  # Gray for outliers
)

# 2. Data Loading and Global Functions
# ------------------------------------------------------------------------------
# - Columns 1-3: Metadata (Column 1: sample ID, Column 2: group1 (CK,Low,High), Column 3: group2 (D0-4))
# - Column 4+: Features (Numeric values for Taxa, Genes, or MGE types)
# - Rows: Each row represents a unique Sample coupled with its Grouping info.
df <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Standard Interquartile Range (IQR) outlier detection
detect_outliers <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

# 3. Scenario A: Consolidated Grouping (e.g., D0 vs. D1-D3 Combined)
# ------------------------------------------------------------------------------
analyze_consolidated <- function(target_gene = "tetracenomycin_C") {
  
  df_target <- df %>%
    select(group1, group2, !!sym(target_gene)) %>%
    filter(group1 %in% c("D0", "D1", "D2", "D3")) %>%
    mutate(super_group = ifelse(group1 == "D0", "D0", "Low")) %>%
    mutate(super_group = factor(super_group, levels = c("D0", "Low"))) %>%
    group_by(super_group) %>%
    mutate(
      is_outlier   = detect_outliers(!!sym(target_gene)),
      outlier_lab  = ifelse(is_outlier, as.character(group2), NA),
      point_color  = ifelse(is_outlier, "Outlier", as.character(group1))
    ) %>% ungroup()
  
  # Perform T-test for subtitle annotation
  t_res <- t.test(as.formula(paste(target_gene, "~ super_group")), data = df_target)
  
  p <- ggplot(df_target, aes(x = super_group, y = !!sym(target_gene))) +
    geom_violin(aes(fill = super_group), trim = FALSE, alpha = 0.6, size = 0.3) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = point_color), width = 0.2, size = 2.5, alpha = 0.8) +
    geom_text_repel(aes(label = outlier_lab), na.rm = TRUE, size = 3, fontface = "bold") +
    stat_compare_means(comparisons = list(c("D0", "Low")), method = "t.test", label = "p.signif") +
    scale_fill_manual(values = group2_colors) +
    scale_color_manual(values = group2_colors) +
    theme_classic() +
    labs(title = paste(target_gene, ": Baseline vs. Treatment"),
         subtitle = paste0("T-test: t(", round(t_res$parameter, 1), ") = ", 
                           round(t_res$statistic, 2), ", p = ", format.pval(t_res$p.value, 3)),
         x = "Consolidated Grouping", y = "Abundance", color = "Subgroups")
  
  return(p)
}

# 4. Scenario B: Faceted Gene Comparison (e.g., D0 vs. D4 across Multiple Genes)
# ------------------------------------------------------------------------------
analyze_faceted <- function(target_genes = c("novobiocin", "tetracenomycin_C", "vancomycin"), 
                            groups = c("D0", "D4")) {
  
  df_comp <- df %>%
    select(group1, group2, all_of(target_genes)) %>%
    filter(group1 %in% groups) %>%
    pivot_longer(cols = all_of(target_genes), names_to = "Gene", values_to = "Abundance") %>%
    group_by(Gene, group1) %>%
    mutate(
      is_outlier   = detect_outliers(Abundance),
      outlier_lab  = ifelse(is_outlier, as.character(group2), NA),
      point_color  = ifelse(is_outlier, "Outlier", as.character(group1))
    ) %>% ungroup()
  
  p <- ggplot(df_comp, aes(x = group1, y = Abundance)) +
    geom_violin(aes(fill = group1), trim = FALSE, alpha = 0.6, size = 0.3) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = point_color), width = 0.15, size = 2, alpha = 0.8) +
    geom_text_repel(aes(label = outlier_lab), na.rm = TRUE, size = 2.5, fontface = "bold") +
    facet_wrap(~Gene, scales = "free_y") +
    stat_compare_means(method = "t.test", method.args = list(alternative = "two.sided"), 
                       label = "p.signif", label.x = 1.5) +
    scale_fill_manual(values = group2_colors) +
    scale_color_manual(values = group2_colors) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face = "bold")) +
    labs(title = "Specific Gene Abundance Comparison", x = "Groups", y = "Abundance")
  
  return(p)
}

# 5. Final Output and Console Statistics
# ------------------------------------------------------------------------------
# Generate Plots
p_consolidated <- analyze_consolidated("tetracenomycin_C")
p_faceted      <- analyze_faceted(c("novobiocin", "tetracenomycin_C", "vancomycin"))

# Print to Plots window
print(p_consolidated)
print(p_faceted)

# Run and print console statistics for Scenario B
cat("\n--- Detailed T-test Statistics (Faceted Genes) ---\n")
target_args <- c("novobiocin", "tetracenomycin_C", "vancomycin")
for(gene in target_args){
  res <- t.test(Abundance ~ group1, data = df %>% 
                  select(group1, !!sym(gene)) %>% filter(group1 %in% c("D0", "D4")) %>%
                  pivot_longer(-group1, names_to = "Gene", values_to = "Abundance"))
  cat(paste0(gene, ": t = ", round(res$statistic, 3), ", p = ", format.pval(res$p.value, 4), "\n"))
}