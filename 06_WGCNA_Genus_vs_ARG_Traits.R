# ==============================================================================
# Title: WGCNA for Microbial Genus Abundance & ARG Traits
# Description: Identifying modules of co-occurring genera and correlating them 
#              with ARG7 and ARG5 trait data.
# Parameters: Power=6, SD_Thresh=1e-12, Merge_Thresh=0.25, Min_Module=20
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(WGCNA)
library(tidyverse)

# Enable multi-threading for performance
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# 2. Data Loading
# ------------------------------------------------------------------------------
# Loading Genus Abundance (id as rows, samples as columns)
df_genus <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Loading Phenotypic Data (Sample as rows, Traits as columns)
df_trait <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# 3. Data Transformation & Synchronization
# ------------------------------------------------------------------------------
# Transpose genus data: WGCNA requires samples as rows
data_expr <- df_genus %>%
  column_to_rownames(colnames(df_genus)[1]) %>%
  t() %>%
  as.data.frame()

# Prepare trait data
pheno_data <- df_trait %>%
  column_to_rownames(colnames(df_trait)[1])

# Match samples between expression and trait data
common_samples <- intersect(rownames(data_expr), rownames(pheno_data))
data_expr <- data_expr[common_samples, ]
pheno_data <- pheno_data[common_samples, ]

# Filter features based on standard deviation (Threshold: 1e-12)
gene_sd <- apply(data_expr, 2, sd)
data_expr <- data_expr[, gene_sd > 1e-12]

cat("Preprocessing finished: ", nrow(data_expr), " samples and ", ncol(data_expr), " genera.\n")

# 4. Network Construction (Power = 6)
# ------------------------------------------------------------------------------
softPower <- 6
# Adjacency and Topological Overlap Matrix (TOM)
adjacency <- adjacency(data_expr, power = softPower, type = "signed")
TOM       <- TOMsimilarity(adjacency)
dissTOM   <- 1 - TOM

# 5. Module Identification and Merging
# ------------------------------------------------------------------------------
# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Initial Module detection (MinSize = 20)
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = 20)
dynamicColors <- labels2colors(dynamicMods)

# Merging modules with similar eigengenes (Threshold = 0.25)
MEList <- moduleEigengenes(data_expr, colors = dynamicColors)
MEs    <- MEList$eigengenes
merge  <- mergeCloseModules(data_expr, dynamicColors, cutHeight = 0.25, verbose = 3)

mergedColors <- merge$colors
mergedMEs    <- merge$newMEs
mergedMEs    <- orderMEs(mergedMEs)

# 6. Result Visualization
# ------------------------------------------------------------------------------
# Plot 1: Cluster Dendrogram
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic", "Merged"), 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Clustering of Genera and Module Merging")

# Plot 2: Module-Trait Correlation Heatmap
moduleTraitCor <- cor(mergedMEs, pheno_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(data_expr))

# Creating significance text for the heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Print Heatmap to Plots Panel
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(pheno_data), 
  yLabels = names(mergedMEs),
  ySymbols = names(mergedMEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-Trait Correlation (Signed Network)"
)

# 7. Data Export
# ------------------------------------------------------------------------------
# Exporting results for supplementary tables
final_membership <- data.frame(Genus = colnames(data_expr), Module = mergedColors)
write_csv(final_membership, "WGCNA_Final_Module_Membership.csv")
write_csv(as.data.frame(mergedMEs) %>% rownames_to_column("SampleID"), "WGCNA_Module_Eigengenes.csv")