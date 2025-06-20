# Load required packages
library(WGCNA)   
library(limma)   
library(ggplot2)  

# Enable multi-threading for faster computation
options(stringsAsFactors = FALSE)
enableWGCNAThreads()


# Input files:
data_expr <- read.table(file.choose(), header = TRUE, sep = "\t", row.names = 1)
pheno_data <- read.table(file.choose(), header = TRUE, sep = "\t", row.names = 1)


# Filter low-expressed genes (SD threshold = 1e-12)
gene_sd <- apply(data_expr, 1, sd)
keep_genes <- gene_sd > 1e-12
data_expr <- data_expr[keep_genes, ]

# Define soft-thresholding powers to test
powers <- c(1:10, seq(from = 12, to=30, by=2))

# Network topology analysis
sft <- pickSoftThreshold(data_expr, 
                         powerVector = powers, 
                         networkType = "signed", 
                         verbose = 5)
sft$powerEstimate  # Recommended soft-thresholding power

# Plot scale-free topology and mean connectivity
par(mfrow = c(1, 2))  # 1 row, 2 columns

# Left: Scale-free topology fit
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, 
     col = "red")

# Right: Mean connectivity
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     col = "red")

# Select soft-thresholding power
softPower <- 6  # Based on scale-free topology criterion

# Construct adjacency matrix
adjacency <- adjacency(data_expr, power = softPower)

# Calculate topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM  # Dissimilarity matrix

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module identification
minModuleSize <- 20  # Minimum module size
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,  # Module splitting sensitivity (0-4)
  pamRespectsDendro = FALSE,  # Flexible branch cutting
  minClusterSize = minModuleSize
)

# Assign module colors
dynamicColors <- labels2colors(dynamicMods)

# Plot dendrogram with module colors
plotDendroAndColors(
  geneTree,
  dynamicColors,
  groupLabels = "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  main = "Gene clustering and module assignment"
)

# Calculate module eigengenes
MEList <- moduleEigengenes(data_expr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Check for missing values
if(any(is.na(MEs))) {
  stop("Module eigengenes contain NA values - check input data or adjust merge threshold")
}

# Calculate module dissimilarity
MEsDiss <- 1 - cor(MEs)

# Ensure ARG trait is numeric
if(!is.numeric(pheno_data$characteristics_ch1)) {
  pheno_data$characteristics_ch1 <- as.numeric(as.factor(pheno_data$characteristics_ch1))
}

# Module-trait correlation analysis
moduleTraitCor <- cor(MEs, pheno_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = ncol(data_expr))

# Create labeled heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = c("ARG7", "ARG5"),  # Phenotypic traits
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7  # Text size adjustment
)