##########################################################################################
## Create gene x gene co-expression matrix containing proportionality metric (rho)
## using human cortex gene expression matrix from Velmeshev 2019 (scRNA-seq)
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(Matrix)
library(propr)

# Start the clock
ptm <- proc.time()


# Ensembl protein-coding genes
ensemblDf <- read.table('../data/Homo_sapiens.GRCh38.84.GeneAnnotations.txt',
	header=T,sep='\t')
protGenes <- subset(ensemblDf,Biotype=='protein_coding')$GeneName
length(protGenes) # 19826

# row index to gene name mapping table
# NOTE: input file below not included in GitHub
mapDf <- read.table('../data/Velmeshev2019_genes.tsv',header=F,sep='\t')
#rownames(mapDf) # check this is default row indices
colnames(mapDf) <- c('EnsemblID','GeneName')
#dim(mapDf) # 65217 genes

mapDf <- subset(mapDf,GeneName %in% protGenes) # keep protein-coding genes
mapDf <- subset(mapDf,!duplicated(GeneName)) # remove duplicated gene names
dim(mapDf) # 19503 genes


# read in gene x sample sparse matrix
# NOTE: input file below not included in GitHub
mat <- readMM('../data/Velmeshev2019_matrix.mtx.gz')
#mat@Dim # 65217 genes x 104559 samples

# filter out: non-pc genes, duplicated gene names, genes with zero counts in all samples
mat <- mat[as.numeric(rownames(mapDf)),]
dimnames(mat) <- list(mapDf$GeneName,1:ncol(mat))

mat <- mat[rowSums(mat)>0,]
#mat@Dim # 18138 genes x 104559 samples
#str(mat)

# transpose matrix so row = sample and column = gene
mat <- t(mat)
#mat@Dim

# convert to normal matrix, +1 to all counts to avoid taking log of zeros below
mat <- as.matrix(mat) + 1

### test run to estimate time
#ptm <- proc.time()
#testCol <- 1000
#rhoMat <- propr::perb(mat[,1:testCol], select = colnames(mat)[1:testCol])@matrix
###

# calculate proportionality matrix (rho)
rhoMat <- propr::perb(mat, select = colnames(mat))@matrix

# output rho matrix
rhoMat <- signif(rhoMat,digits=3)
write.table(rhoMat,gzfile('../output/Velmeshev2019_CoExpPropRho.txt.gz'),sep='\t',quote=F)

# Stop the clock
proc.time() - ptm

