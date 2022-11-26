##########################################################################################
## Generate gene expression heat map for:
## (1) 24 ASD risk genes from Satterstrom et al.
## (2) housekeeping or cell type marker genes
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

#BiocManager::install("edgeR")
library(edgeR)
library(reshape2)
library(ggplot2)


# ----------------------------------------------------------------------------------------
# read in gene count matrix derived from bulk RNA-seq data
counts <- read.table('../data/ASD_RNAseq_geneCounts_noTail.tsv',
	header=T,sep='\t',row.names=1)
counts <- counts[rowSums(counts==0)<6,] # only store genes with > 1 non-zero counts
dim(counts) # 37857 genes (w/ >1 non-zero counts) x 6 samples
sampLabels <- c('d0-r1','d0-r2','d21-r1','d21-r2','d51-r1','d51-r2')
names(counts) <- sampLabels

# create  DGEList object
dgeObj <- DGEList(counts)

# sample annotations
group <- as.factor(rep(c('d0','d21','d51'), c(2,2,2)))
dgeObj$samples$group <- group
dgeObj$samples

# ENSG ID to HGNC symbol mapping
genes <- read.table('../data/Homo_sapiens.GRCh38.84.GeneAnnotations.txt',
	header=T,sep='\t',stringsAsFactors=F)
dim(genes) # same # of genes as in count matrix

# arrange genes in same order as count matrix
genes <- genes[match(rownames(counts),genes$GeneID),]
dgeObj$genes <- genes


# filter low-count genes 
keep.exprs <- filterByExpr(dgeObj,min.count=5) 
dgeObj <- dgeObj[keep.exprs,, keep.lib.sizes=FALSE]
dim(dgeObj) # 22269 genes (14845 protein_coding) 


# calculate normalized gene expression using TMM (trimmed mean of M-vlaues)
dgeObj <- calcNormFactors(dgeObj, method='TMM')
#samples$norm.factors

# get TMM-normalized values
lcpm_norm <- cpm(dgeObj,log=T)
rownames(lcpm_norm) <- dgeObj$genes$GeneName


# calculate avg. log2 CPM for each time point
lcpm_avg <- data.frame(Gene=rownames(lcpm_norm))
for (i in c(0,21,51)) {
	r1 <- paste('d',i,'-r1',sep='')
	r2 <- paste('d',i,'-r2',sep='')

	avgCol <- paste('d',i,sep='')
	lcpm_avg[,avgCol] <- (lcpm_norm[,r1] + lcpm_norm[,r2])/2
}

# output avg. gene expression matrix
write.table(lcpm_avg,'../output/ASD_RNA_geneAvgLcpm.txt',sep='\t',quote=F,row.names=F)


# -------------------------------------------------------------------------------------------
# heat map of 24 ASD Satterstrom genes
x <- c('ADNP','ANK2','ANKRD11','ARID1B','ASH1L','CHD2','CHD8','CTNNB1',
	'DEAF1','DSCAM','DYRK1A','FOXP1','GIGYF1','GRIN2B','KDM6B','KMT5B',
	'MED13L','POGZ','PTEN','SCN2A','SHANK3','SLC6A1','SYNGAP1','TLK2')

df <- lcpm_avg[match(x,lcpm_avg$Gene),]
colnames(df) <- c('Gene','week 0','week 3','week 6')

# data frame for plotting
meltDf <- melt(df,id='Gene')
meltDf$Gene <- factor(meltDf$Gene,levels=rev(x))
meltDf$variable <- factor(meltDf$variable,levels=unique(meltDf$variable))

# cap max value at 9 for better color spread
meltDf$value[meltDf$value > 9] <- 9

# heat map
pdf('../output/ASD_RNA_SatterstromHeatMap.pdf',height=3.5,width=2.5)

ggplot(meltDf,aes(x=variable,y=Gene,fill=value)) + geom_tile() +

scale_fill_distiller(name=bquote(log[2]*'(CPM)'),palette='RdYlBu',
	direction=-1,breaks=c(0,4,8)) +

xlab('Time point') + ylab('Gene') + theme_minimal() +
theme(legend.position='bottom',legend.title=element_text(size=8),
	legend.text=element_text(size=8),legend.key.size=unit(0.5,'line'),
	legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5),
	axis.title=element_text(size=8),axis.text=element_text(size=8))

dev.off()


# -------------------------------------------------------------------------------------------
# heat map of housekeeping or marker genes
markerDf <- read.table('../data/Nehme2018_MarkerGenes.txt',header=T,sep='\t',stringsAsFactors=F)
markerDf <- subset(markerDf,MarkerType!='Progenitor')
markerDf$Color <- 'black'
markerDf$Color[markerDf$MarkerType=='Pluripotency'] <- 'darkcyan'
markerDf$Color[markerDf$MarkerType=='Neuronal'] <- 'purple'

df <- lcpm_avg[match(markerDf$Gene,lcpm_avg$Gene),]
colnames(df) <- c('Gene','week 0','week 3','week 6')

# data frame for plotting
meltDf <- melt(df,id='Gene')
meltDf$Gene <- factor(meltDf$Gene,levels=rev(df$Gene))
meltDf$variable <- factor(meltDf$variable,levels=unique(meltDf$variable))

# heat map
pdf('../output/ASD_RNA_MarkerGenesHeatMap.pdf',height=2.5,width=2.2)

ggplot(meltDf,aes(x=variable,y=Gene,fill=value)) + geom_tile() +
scale_fill_distiller(name=bquote(log[2]*'(CPM)'),palette='RdYlBu',direction=-1) +
xlab('Time point') + ylab('Gene') + theme_minimal() +
theme(legend.position='bottom',legend.title=element_text(size=8),legend.text=element_text(size=8),
	legend.key.size=unit(0.5,'line'),
	legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5),
	axis.title=element_text(size=8),axis.text.x=element_text(size=8),
	axis.text.y=element_text(size=8,color=rev(markerDf$Color)))

dev.off()


