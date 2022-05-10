##########################################################################################
## Create gene x gene co-expression matrix containing proportionality metric (rho)
## using BrainSpan RNA-seq data
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

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
rowDf <- read.csv('../data/BrainSpan_rows_metadata.csv',header=T,row.names=1)
#rownames(rowDf) # check this is default row indices
#dim(rowDf) # 52376 genes

rowDf <- subset(rowDf,gene_symbol %in% protGenes) # keep protein-coding genes
rowDf <- subset(rowDf,!duplicated(gene_symbol)) # remove duplicated gene names
dim(rowDf) # 18090 genes


# read in gene x sample matrix
# NOTE: input file below not included in GitHub
mat <- read.csv('../data/BrainSpan_expression_matrix.csv',header=F,row.names=1)
colnames(mat) <- 1:ncol(mat)
#dim(mat) # 52375 genes

# filter out: non-pc genes, duplicated gene names, genes with zero counts in all samples
mat <- mat[as.numeric(rownames(rowDf)),]
rownames(mat) <- rowDf$gene_symbol

mat <- mat[rowSums(mat)>0,]
#dim(mat) # 17975 genes 

# transpose matrix so row = sample and column = gene
mat <- t(mat)
dim(mat)

# +1 to all counts to avoid taking log of zeros below
mat <- mat + 1

### test run to estimate time
#ptm <- proc.time()
#testCol <- 1000
#rhoMat <- propr::perb(mat[,1:testCol], select = colnames(mat)[1:testCol])@matrix
###

# calculate proportionality matrix (rho)
rhoMat <- propr::perb(mat, select = colnames(mat))@matrix

# output rho matrix
rhoMat <- signif(rhoMat,digits=3)
write.table(rhoMat,gzfile('../output/BrainSpan_CoExpPropRho.txt.gz'),sep='\t',quote=F)

# Stop the clock
proc.time() - ptm

