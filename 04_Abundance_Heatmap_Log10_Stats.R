# ==============================================================================
# Title: Feature Abundance Heatmap
# Description: Generates heatmaps and exports a complete 
#              matrix including Mean, Log10-values, and P-values for all groups.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(pheatmap)

# 2. Configuration & Offsets
# ------------------------------------------------------------------------------
log_offset <- 1e-6  # Transformation: log10(Abundance + 1e-6)
feature_label <- "ARG" # Label prefix for the output

# 3. Data Loading and Preprocessing
# ------------------------------------------------------------------------------
# - Columns 1-4: Metadata (Column 1: sample ID, Column 2: group1 (CK,Low,High), Column 3: group2 (D0-4), Column 4: Mixture (specific pesticide combination)
# - Column 5+: Features (Numeric values for Taxa, Genes, or MGE types)
# - Rows: Each row represents a unique Sample coupled with its Grouping info.
df <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
val_cols <- df %>% select(where(is.numeric)) %>% colnames()

df_long <- df %>%
  pivot_longer(cols = all_of(val_cols), names_to = "FeatureClass", values_to = "Abundance")

# Set consistent group order
all_groups <- unique(df$group1)
target_order <- c("CK", intersect(c("Low", "High"), all_groups), setdiff(all_groups, c("CK", "Low", "High")))

# 4. Systematic Abundance & Statistical Calculation
# ------------------------------------------------------------------------------
# A. Calculate Mean Abundance for ALL groups
summary_table <- df_long %>%
  group_by(FeatureClass, group1) %>%
  summarise(mean_abd = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  # Add log-transformed values for all groups
  mutate(log_mean = log10(mean_abd + log_offset))

# B. Calculate P-values (Treatment vs. CK)
# We perform this separately to merge into the final wide table
test_groups <- setdiff(all_groups, "CK")
p_values <- data.frame()

for (g in test_groups) {
  for (f in unique(df_long$FeatureClass)) {
    val_test <- df_long %>% filter(FeatureClass == f, group1 == g) %>% pull(Abundance)
    val_ref  <- df_long %>% filter(FeatureClass == f, group1 == "CK") %>% pull(Abundance)
    
    p <- tryCatch({
      if(sd(val_test) == 0 && sd(val_ref) == 0) 1
      else t.test(val_test, val_ref)$p.value
    }, error = function(e) 1)
    
    p_values <- rbind(p_values, data.frame(FeatureClass = f, group1 = g, p_val = p))
  }
}

# 5. Integrate into Final Wide Table for Export
# ------------------------------------------------------------------------------
# Create a professional wide table containing CK, Low, High data in one row per feature
export_table <- summary_table %>%
  pivot_wider(names_from = group1, values_from = c(mean_abd, log_mean)) %>%
  left_join(
    p_values %>% pivot_wider(names_from = group1, values_from = p_val, names_prefix = "p_val_vs_CK_"),
    by = "FeatureClass"
  )

# Export the comprehensive table
write.csv(export_table, paste0(feature_label, "_Full_Heatmap_Data_Report.csv"), row.names = FALSE)

# 6. Prepare Heatmap Matrix and Significance Stars
# ------------------------------------------------------------------------------
# Matrix for colors (Mean values)
heatmap_mat <- export_table %>%
  select(FeatureClass, starts_with("mean_abd_")) %>%
  column_to_rownames("FeatureClass")
colnames(heatmap_mat) <- gsub("mean_abd_", "", colnames(heatmap_mat))
heatmap_mat <- heatmap_mat[, target_order, drop = FALSE]

# Matrix for stars (P-values)
sig_labels <- heatmap_mat
sig_labels[] <- ""

get_stars <- function(p) {
  if (is.na(p) | p >= 0.05) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  return("*")
}

# Fill stars for treatment groups (CK remains empty)
for (g in test_groups) {
  p_col <- paste0("p_val_vs_CK_", g)
  for (f in rownames(sig_labels)) {
    p_val <- export_table[[p_col]][export_table$FeatureClass == f]
    sig_labels[f, g] <- get_stars(p_val)
  }
}

# 7. Figure Generation
# ------------------------------------------------------------------------------
plot_mat <- log10(heatmap_mat + log_offset)

final_heatmap <- pheatmap(
  plot_mat,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = sig_labels, 
  fontsize_number = 12,
  number_color = "black",
  color = colorRampPalette(c("#f7fbff", "#9ecae1", "#084594"))(100),
  border_color = "white",
  main = paste("Mean", feature_label, "Abundance (log10 Scale)"),
  angle_col = 0
)