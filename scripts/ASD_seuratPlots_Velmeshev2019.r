##########################################################################################
## Generate seurat tSNE plots for ASD baits using scRNA-seq data from Velmeshev2019
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)

# ----------------------------------------------------------------------------------------
# read in Velmeshev2019 data

# read in normalized, log2-transformed expression matrix
# NOTE: input file below not included in GitHub
expMatrix <- fread('zcat < ../data/Velmeshev2019_exprMatrix.tsv.gz')
genes <- expMatrix[,1][[1]]
expMatrix = data.frame(expMatrix[,-1], row.names=genes)
dim(expMatrix) # 36501 genes x 104559 cells

# read in meta data (including cell type cluster info)
# NOTE: input file below not included in GitHub
meta <- read.table('../data/Velmeshev2019_meta.tsv',header=T,sep='\t',as.is=T,row.names=1)
rownames(meta) <-  gsub('-','.',rownames(meta))

# craete seurat object
seuratObj <- CreateSeuratObject(counts=expMatrix, meta.data=meta)

# set cluster ID for each cell
Idents(seuratObj) <- 'cluster'
Idents(seuratObj) <- factor(Idents(seuratObj),levels=sort(levels(Idents(seuratObj))))

# set tSNE coordinates for each cell
# NOTE: input file below not included in GitHub
tsneMatrix <- data.matrix(read.table('../data/Velmeshev2019_tSNE.coords.tsv.gz',
	header=F,sep='\t',row.names=1))
colnames(tsneMatrix) <- c('tSNE_1','tSNE_2')
rownames(tsneMatrix) <- gsub('-','.',rownames(tsneMatrix))
seuratObj[['tsne']] <- CreateDimReducObject(tsneMatrix,key='tSNE_',assay='RNA')


# ----------------------------------------------------------------------------------------
# tSNE plots

# tSNE plot with cluster name labels
png("../output/ASD_Velmeshev2019_DimPlot.png", width=5, height=5, units="in", res=300)

DimPlot(seuratObj,label=T,label.size=2.5) +
theme(axis.title=element_text(size=8),axis.text=element_blank(),
	axis.ticks=element_blank(),legend.position='none')

dev.off()


# tSNE plots for each of the 13 baits
baitTable <- data.frame(
	Gene=c('ADNP','ANK2','ARID1B','CHD8','CTNNB1','DYRK1A','GIGYF1',
		'MED13L','PTEN','SCN2A','SHANK3','SYNGAP1','TLK2'),
	ID=c('ENSG00000101126','ENSG00000145362','ENSG00000049618','ENSG00000100888',
		'ENSG00000168036','ENSG00000157540','ENSG00000146830','ENSG00000123066',
		'ENSG00000171862','ENSG00000136531','ENSG00000251322','ENSG00000197283',
		'ENSG00000146872'))

png('../output/ASD_Velmeshev2019_BaitFeaturePlots.png',
	height=6,width=11,units="in",res=300)
	
fPlots <- FeaturePlot(seuratObj, features=as.character(baitTable$ID), 
	cols=c('lightgrey','red'), combine=F)

for (i in 1:length(fPlots)) {
	fPlots[[i]] <- fPlots[[i]] + ggtitle(baitTable$Gene[baitTable$ID==baitTable$ID[i]]) +
		theme(plot.title=element_text(size=10,face='bold'),
			axis.title=element_text(size=8),axis.text=element_blank(),
			axis.ticks=element_blank(),
			legend.text=element_text(size=8),legend.key.width=unit(0.6,"line"))
}

CombinePlots(fPlots,ncol=5)

dev.off()

