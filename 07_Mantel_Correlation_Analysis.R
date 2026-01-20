# ==============================================================================
# Title: Mantel Test and Correlation Matrix Visualization
# Description: Mantel test links and Pearson correlation heatmap.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(linkET)
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)

# 2. Data Loading & Matrix Preprocessing
# ------------------------------------------------------------------------------
# Loading ARG data (Environment/Target profiles)
# Columns: sample, specific ARG abundance
df_arg_raw <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
df_arg <- df_arg_raw %>% column_to_rownames(colnames(df_arg_raw)[1])

# Loading MGE data (Predictor/Source profiles)
# Columns: sample, specific MGE abundance
df_mge_raw <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
df_mge <- df_mge_raw %>% column_to_rownames(colnames(df_mge_raw)[1])

# Synchronize Samples between datasets
common_samples <- intersect(rownames(df_arg), rownames(df_mge))
df_arg <- df_arg[common_samples, ]
df_mge <- df_mge[common_samples, ]

# Automatically adapt to all columns in the MGE table
all_mges <- colnames(df_mge)
spec_list <- setNames(as.list(1:length(all_mges)), all_mges)

# 3. Mantel Test Calculation & Significance Mapping
# ------------------------------------------------------------------------------
mantel_result <- mantel_test(df_mge, df_arg, spec_select = spec_list) %>% 
  mutate(
    # Mapping line width to Mantel's r (Effect Size)
    rd = cut(abs(r), breaks = c(-Inf, 0.1, 0.3, Inf),
             labels = c("< 0.1", "0.1 - 0.3", ">= 0.3")),
    
    # Mapping line type to P-value (Significance)
    pd = cut(p, breaks = c(-Inf, 0.05, Inf),
             labels = c("< 0.05", ">= 0.05")),
    
    # Mapping line color to Correlation Direction (Positive/Negative)
    r_sign = ifelse(r > 0, "Positive", "Negative")
  )

# 4. Core Visualization (Publication Quality)
# ------------------------------------------------------------------------------
# Calculate Pearson correlation for the background matrix
cor_matrix <- cor(df_arg)

# Define professional color palettes based on reference image
# Pearson Heatmap: Deep Blue -> Light Yellow -> Deep Red
pearson_pal <- c("#445494", "#7997C0", "#F7E092", "#DD6D46", "#B2182B")
# Mantel Links: Positive = Red, Negative = Blue
link_pal    <- c("Positive" = "#B2182B", "Negative" = "#5681B9")

# Constructing the plot
p <- qcorrplot(cor_matrix, 
               type = "upper", 
               diag = FALSE, 
               grid.col = "grey95", 
               grid.size = 0.1) +
  # Layer 1: Background Correlation Square
  geom_square() +
  
  # Layer 2: Foreground Mantel Test Links
  geom_couple(aes(colour = r_sign, size = rd, linetype = pd), 
              data = mantel_result, 
              curvature = 0.15) +
  
  # Color Scale for Pearson Matrix
  scale_fill_gradientn(colours = pearson_pal,
                       limits = c(-0.2, 0.4), 
                       name = "Pearson's r") +
  
  # Color Scale for Mantel's Direction
  scale_colour_manual(values = link_pal, 
                      name = "Mantel's r sign") +
  
  # Size Scale for Mantel's Strength
  scale_size_manual(values = c("< 0.1" = 0.4, "0.1 - 0.3" = 1.2, ">= 0.3" = 2.5), 
                    name = "Mantel's r abs") +
  
  # Line Type Scale for Significance
  scale_linetype_manual(values = c("< 0.05" = "solid", ">= 0.05" = "dashed"), 
                        name = "Mantel's P") +
  
  # Final Polishing and Legend Optimization
  guides(
    fill = guide_colorbar(order = 1, barwidth = 0.8, barheight = 5),
    colour = guide_legend(order = 2, override.aes = list(linewidth = 1.5)),
    linetype = guide_legend(order = 3),
    size = guide_legend(order = 4)
  ) +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(face = "bold", size = 9),
    legend.key = element_blank(),
    panel.background = element_rect(fill = "white")
  )

# 5. Output and Export
# ------------------------------------------------------------------------------
# Display plot in RStudio
print(p)

# Export
write_csv(mantel_result, "Mantel_Statistical_Summary.csv")