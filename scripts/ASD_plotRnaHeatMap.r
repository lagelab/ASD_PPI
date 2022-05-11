##########################################################################################
## Generate gene expression heat map for 24 ASD risk genes from Satterstrom et al.
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
counts <- counts[rowSums(counts==0)<6,] # only store genes with >1 non-zero counts
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


# calculate normalization factors using TMM method (trimmed mean of M-vlaues)
dgeObj <- calcNormFactors(dgeObj, method='TMM')
#samples$norm.factors

# get TMM-normalized values
lcpm_norm <- cpm(dgeObj,log=T)
rownames(lcpm_norm) <- dgeObj$genes$GeneName


# ----------------------------------------------------------------------------------------
# create heat map of 24 Satterstrom genes

x <- c('ADNP','ANK2','ANKRD11','ARID1B','ASH1L','CHD2','CHD8','CTNNB1',
	'DEAF1','DSCAM','DYRK1A','FOXP1','GIGYF1','GRIN2B','KDM6B','KMT5B',
	'MED13L','POGZ','PTEN','SCN2A','SHANK3','SLC6A1','SYNGAP1','TLK2')

df <- as.data.frame(lcpm_norm[match(x,rownames(lcpm_norm)),])

# calculate avg. log2 CPM for each time point
for (i in c(0,21,51)) {
	r1 <- paste('d',i,'-r1',sep='')
	r2 <- paste('d',i,'-r2',sep='')

	avgCol <- paste('d',i,sep='')
	df[,avgCol] <- (df[,r1]+df[,r2])/2
}
df$Gene <- rownames(df)
df <- df[,c('Gene','d0','d21','d51')]
colnames(df) <- c('Gene','week 0','week 3','week 6')

# data frame for plotting
meltDf <- melt(df,id='Gene')
meltDf$Gene <- factor(meltDf$Gene,levels=rev(x))
meltDf$variable <- factor(meltDf$variable,levels=unique(meltDf$variable))

# cap max value at 9 for better color spread
meltDf$value[meltDf$value > 9] <- 9

# heat map
pdf('../output/ASD_RnaHeatMap.pdf',height=3.5,width=2.5)

ggplot(meltDf,aes(x=variable,y=Gene,fill=value)) + geom_tile() +

scale_fill_distiller(name=bquote(log[2]*'(CPM)'),palette='RdYlBu',
	direction=-1,breaks=c(0,4,8)) +

xlab('Time point') + ylab('Gene') + theme_minimal() +
theme(legend.position='bottom',legend.title=element_text(size=8),
	legend.text=element_text(size=8),legend.key.size=unit(0.5,'line'),
	legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5),
	axis.title=element_text(size=8),axis.text=element_text(size=8))

dev.off()

