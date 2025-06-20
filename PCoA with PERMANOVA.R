# Load required libraries
library(vegan) 

# Import abundance data
df <- read.table(file.choose(), header = TRUE, sep = '\t', row.names = 1)

# Import metadata with grouping information
# First column should match sample names in abundance data
data1 <- read.table(file.choose(), sep = '\t', header = TRUE, row.names = 1)

# Calculate Bray-Curtis dissimilarity matrix
distance <- vegdist(df, method = 'bray')

# Perform Principal Coordinates Analysis (PCoA)
# k = number of dimensions (here set to n-1 where n is number of samples)
# eig = TRUE returns eigenvalues
pcoa <- cmdscale(distance, k = (nrow(df) - 1), eig = TRUE)

# Convert distance matrix to regular matrix format
distance_matrix <- as.matrix(distance)

# Perform PERMANOVA (Permutational Multivariate Analysis of Variance)
# Tests for significant differences between groups
# permutations = 999: number of permutations for significance testing
permanova_result <- adonis2(df ~ data1$group, 
                            data = df,
                            permutations = 999, 
                            method = 'bray')

# Display PERMANOVA results
print(permanova_result)
