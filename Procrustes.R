
# Libraries
library(vegan)

# Load data files
df<- read.table(file.choose(),header=T,sep='\t',row.names=1)
env <- read.table(file.choose(),sep="\t",header = T,row.names = 1,check.names = F)

# Calculate dissimilarity matrices
spe.dist <- vegdist(df, method = "bray")
env.dist <- vegdist(scale(env), "euclid")

# Ordination analysis
set.seed(123)
mds.s <- monoMDS(spe.dist)
mds.e <- monoMDS(env.dist)

# Procrustes analysis
pro.s.e <- procrustes(mds.s,mds.e, symmetric = TRUE)
summary(pro.s.e)