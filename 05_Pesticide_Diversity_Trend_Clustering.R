# ==============================================================================
# Title: Pesticide Diversity Gradient Clustering (Mfuzz)
# Description: Fuzzy c-means clustering to identify response patterns of 
#              pesticide diversity across a treatment gradient (D0-D4).
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(Mfuzz)
library(tidyverse)

# 2. Data Loading and Preprocessing
# ------------------------------------------------------------------------------
# Input: Rows are pesticide features, Columns are gradient points (D0, D1, D2, D3, D4)
# Note: Data should be mean values per gradient point.
raw_data <- read.table(file.choose(), header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
pesticide_matrix <- as.matrix(raw_data)

# Create ExpressionSet object
mfuzz_set <- new('ExpressionSet', exprs = pesticide_matrix)

# 3. Data Refinement & Standardization
# ------------------------------------------------------------------------------
# A. Filter features with excessive missing values or zero variance
mfuzz_set <- filter.NA(mfuzz_set, thres = 0.25)
mfuzz_set <- fill.NA(mfuzz_set, mode = "mean")
mfuzz_set <- filter.std(mfuzz_set, min.std = 0)

# B. Z-score Standardization
# Essential for grouping pesticides with similar response 'shapes' regardless of absolute scale
mfuzz_set_std <- standardise(mfuzz_set)

# 4. Fuzzy C-Means Clustering Logic
# ------------------------------------------------------------------------------
# A. Estimate the fuzzifier parameter 'm'
m_value <- mestimate(mfuzz_set_std)

# B. Number of Clusters
set.seed(123) 
cluster_num <- 8

# C. Execute clustering
mfuzz_res <- mfuzz(mfuzz_set_std, c = cluster_num, m = m_value)

# 5. Visualization
# ------------------------------------------------------------------------------
# Adjusting graphical parameters for manuscript submission

mfuzz.plot2(
  mfuzz_set_std, 
  cl = mfuzz_res, 
  mfrow = c(2, 4),           # Layout for 8 clusters
  time.labels = colnames(pesticide_matrix), 
  xlab = "Diversity Gradient", 
  ylab = "Standardized Diversity Index",
  col.lab = "black", 
  col.main = "darkred",      # Distinction for pesticide diversity
  family = "sans", 
  cex.main = 1.3,
  x11 = FALSE
)

dev.off()

# 6. Data Export: Comprehensive Report
# ------------------------------------------------------------------------------
# Extract cluster information and membership scores
cluster_info <- data.frame(
  Pesticide_ID = names(mfuzz_res$cluster),
  Cluster = mfuzz_res$cluster,
  Membership_Score = apply(mfuzz_res$membership, 1, max)
)

# Merge with original data for a complete summary
final_report <- merge(raw_data, cluster_info, by.x = "row.names", by.y = "Pesticide_ID")
colnames(final_report)[1] <- "Pesticide_Feature"

# Export the detailed report
write.csv(final_report, "Pesticide_Gradient_Cluster_Report.csv", row.names = FALSE)