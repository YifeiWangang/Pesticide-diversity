# Install Mfuzz package from Bioconductor
# install.packages('BiocManager')
# BiocManager::install('Mfuzz')

# Load Mfuzz package
library(Mfuzz)

# Read expression matrix file (biological replicates should be averaged beforehand)
gene_data <- read.table(file.choose(), header=TRUE, sep='\t', row.names=1)
gene_data <- as.matrix(gene_data)

# Create ExpressionSet object
mfuzz_class <- new('ExpressionSet', exprs = gene_data)

# Standardize data (z-score normalization)
mfuzz_class <- standardise(mfuzz_class)

# Cluster analysis using fuzzy c-means algorithm (see ?mfuzz for details)
# User must define the number of clusters (here set to 8)
# Set random seed for reproducibility
set.seed(123)
cluster_num <- 8
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

# Plot cluster profiles (see ?mfuzz.plot2 for details)
# time.labels parameter sets sample/time point labels
# Additional parameters can adjust colors, line width, axes, fonts, etc.
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 4), 
            time.labels = colnames(gene_data))

# Check number of genes in each cluster
cluster_size <- mfuzz_cluster$size
names(cluster_size) <- 1:cluster_num
cluster_size

# View cluster assignment for each gene (showing first 6 entries)
head(mfuzz_cluster$cluster)