# ==============================================================================
# Title: Beta-diversity Analysis (PCoA) and Statistical Distance Testing (T-test)
# Description: Generate PCoA coordinates, centroids, and calculate statistical 
#          significance of shifts from baseline (D0) using t-tests.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggpubr)

# Define color palette for overall metagenome
group2_colors <- c(
  "D0" = "#7B81AF",
  "D1" = "#D3D3D3",
  "D2" = "#73AB81",
  "D3" = "#87619D",
  "D4" = "#72A6A5"
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
# - Columns 1-3: Metadata (Column 1: sample ID, Column 2: group1 (CK,Low,High), Column 3: group2 (D0-4))
# - Column 4+: Features (Numeric values for Taxa, Genes, or MGE types)
# - Rows: Each row represents a unique Sample coupled with its Grouping info.
df <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Construct community matrix
data_matrix <- df %>% 
  column_to_rownames("sample") %>% 
  select(-group1, -group2)

# Compute Bray-Curtis dissimilarity matrix
dist_matrix <- vegdist(data_matrix, method = "bray")

# 3. PCoA Calculation and Export of Coordinates/Centroids
# ------------------------------------------------------------------------------
# Perform Principal Coordinate Analysis
pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# Assemble PCoA data frame and merge with metadata
pcoa_df <- as.data.frame(pcoa_res$points) %>%
  rename(PCoA1 = V1, PCoA2 = V2) %>%
  rownames_to_column("sample") %>%
  left_join(df %>% select(sample, group1, group2), by = "sample")

# Calculate group centroids for spatial orientation and visualization
centroids <- pcoa_df %>%
  group_by(group2) %>%
  summarise(cX = mean(PCoA1), cY = mean(PCoA2), .groups = "drop")

# Export quantitative data
write.csv(pcoa_df, "PCoA_Coordinates_Export.csv", row.names = FALSE)
write.csv(centroids, "PCoA_Group_Centroids.csv", row.names = FALSE)

# 4. Statistical Analysis: Mean Distance to Control (D0/CK)
# ------------------------------------------------------------------------------
# Isolate Control (D0) samples to calculate distances of treatment groups to baseline
d0_samples <- df$sample[df$group2 == "D0"]
dist_mat_full <- as.matrix(dist_matrix)

# Compute average dissimilarity of each sample to the D0 group replicates
d0_dist_data <- data.frame(
  sample = colnames(dist_mat_full),
  group2 = df$group2,
  dist_to_CK = colMeans(dist_mat_full[d0_samples, , drop = FALSE])
) %>% filter(group2 != "D0")

# Perform Student's t-test for pairwise group comparisons (Unadjusted P-values)
pairwise_t_test <- compare_means(dist_to_CK ~ group2, data = d0_dist_data, 
                                 method = "t.test", p.adjust.method = "none")

# Export statistical test results
write.csv(as.data.frame(pairwise_t_test), "Distance_to_CK_T_Test_Results.csv", row.names = FALSE)

# Global community structure significance (PERMANOVA)
set.seed(123)
adonis_res <- adonis2(dist_matrix ~ group2, data = df)

# 5. Figure 1: PCoA Visualization
# ------------------------------------------------------------------------------
pcoa_df_plot <- pcoa_df %>% left_join(centroids, by = "group2")

pcoa_plot <- ggplot(pcoa_df_plot, aes(x = PCoA1, y = PCoA2, color = group2)) +
  geom_segment(aes(xend = cX, yend = cY), alpha = 0.2, show.legend = FALSE) +
  stat_ellipse(aes(fill = group2), geom = "polygon", alpha = 0.1, linetype = "dashed", level = 0.95) +
  geom_point(size = 3, alpha = 0.8) +
  geom_point(data = centroids, aes(x = cX, y = cY, fill = group2), 
             size = 5, shape = 21, color = "black", stroke = 1.2, show.legend = FALSE) +
  scale_color_manual(values = group2_colors) +
  scale_fill_manual(values = group2_colors) +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = "right") +
  labs(
    x = paste0("PCoA1 (", round(pcoa_res$eig[1]/sum(pcoa_res$eig)*100, 2), "%)"),
    y = paste0("PCoA2 (", round(pcoa_res$eig[2]/sum(pcoa_res$eig)*100, 2), "%)"),
    title = "PCoA of Community Composition"
  ) +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("PERMANOVA: P = ", adonis_res$`Pr(>F)`[1]),
           hjust = 1.1, vjust = -1, size = 3.5, fontface = "italic")

# 6. Figure 2: Distance to Baseline Boxplot
# ------------------------------------------------------------------------------
# Define specific comparisons for t-test annotation
target_comparisons <- list( c("D1", "D4"), c("D2", "D4"), c("D3", "D4") )

dist_box_plot <- ggplot(d0_dist_data, aes(x = group2, y = dist_to_CK, fill = group2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
  # Apply pairwise t-test significance labels (p.signif showing *, **, etc.)
  stat_compare_means(comparisons = target_comparisons, method = "t.test", label = "p.signif") + 
  scale_fill_manual(values = group2_colors) +
  theme_classic() +
  labs(
    x = "Group",
    y = "Bray-Curtis Distance to CK (D0)",
    title = "Dissimilarity to Baseline"
  ) +
  theme(legend.position = "none")

# 7. Final Output
# ------------------------------------------------------------------------------
# Plot Figures individually
print(pcoa_plot)
print(dist_box_plot)