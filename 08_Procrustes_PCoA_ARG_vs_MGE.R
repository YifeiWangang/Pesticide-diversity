# ==============================================================================
# Title: Procrustes Analysis based on PCoA (Principal Coordinate Analysis)
# Description: Using PCoA coordinates to evaluate community congruence.
# ==============================================================================

library(vegan)
library(tidyverse)
library(ggplot2)

# 1. Data Preparation
# ------------------------------------------------------------------------------
# - Rows : Samples
# - Columns : Features (e.g., Antibiotic classes or MGE types)

df_arg <- read_delim(file.choose(), delim = "\t", show_col_types = FALSE) %>% 
  column_to_rownames("sample")

df_mge <- read_delim(file.choose(), delim = "\t", show_col_types = FALSE) %>% 
  column_to_rownames("sample")

common <- intersect(rownames(df_arg), rownames(df_mge))
df_arg <- df_arg[common, ]
df_mge <- df_mge[common, ]

# 2. PCoA Ordination (Using cmdscale)
# ------------------------------------------------------------------------------
# Calculate Bray-Curtis distance
dist_arg <- vegdist(df_arg, method = "bray")
dist_mge <- vegdist(df_mge, method = "bray")

# Perform PCoA (Principal Coordinate Analysis)
pcoa_arg <- cmdscale(dist_arg, k = 2, eig = TRUE)
pcoa_mge <- cmdscale(dist_mge, k = 2, eig = TRUE)

# Extract PCoA coordinates
points_arg <- pcoa_arg$points
points_mge <- pcoa_mge$points

# 3. Procrustes Analysis & Statistical Test
# ------------------------------------------------------------------------------
# Procrustes rotation based on PCoA axes
pro_fit <- procrustes(points_arg, points_mge, symmetric = TRUE)

# Protest for significance
set.seed(123)
pro_test <- protest(points_arg, points_mge, permutations = 999)

# Get precise metrics
m2_actual <- round(pro_test$ss, 3) 
p_actual  <- pro_test$signif

# 4. Visualization
# ------------------------------------------------------------------------------
df_points_arg <- as.data.frame(pro_fit$X)
df_points_mge <- as.data.frame(pro_fit$Yrot)
colnames(df_points_arg) <- colnames(df_points_mge) <- c("PCoA1", "PCoA2")

plot_df <- data.frame(
  x_arg = df_points_arg$PCoA1, y_arg = df_points_arg$PCoA2,
  x_mge = df_points_mge$PCoA1, y_mge = df_points_mge$PCoA2
)

# 
ggplot(plot_df) +
  geom_segment(aes(x = x_arg, y = y_arg, xend = x_mge, yend = y_mge),
               arrow = arrow(length = unit(0.1, "cm")), 
               color = "grey88", size = 0.35) +
  # ARG: Filled Red Circles
  geom_point(aes(x = x_arg, y = y_arg, fill = "ARG"), 
             color = "black", size = 3, shape = 21, stroke = 0.3) +
  # MGE: Open Blue Circles
  geom_point(aes(x = x_mge, y = y_mge, color = "MGE"), 
             size = 3, shape = 1, stroke = 1) +
  scale_fill_manual(values = c("ARG" = "#B2182B"), name = NULL) +
  scale_color_manual(values = c("MGE" = "#445494"), name = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.88, 0.12),
    legend.background = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Procrustes Analysis (PCoA based)",
    x = "PCoA 1", y = "PCoA 2"
  ) +
  # Annotation using M2 and P-value
  annotate("text", x = -Inf, y = Inf, 
           label = paste0("italic(M)^2 == ", m2_actual, "\nitalic(P) == ", p_actual),
           hjust = -0.2, vjust = 1.5, size = 4.5, fontface = "bold", parse = TRUE)