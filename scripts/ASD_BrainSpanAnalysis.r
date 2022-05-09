##########################################################################################
## BrainSpan expression profiles for various ASD gene sets
##
## Author: Kalliopi Tsafou, Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(dplyr)
library(reshape2)
library(ggplot2)


# ----------------------------------------------------------------------------------------
# BrainSpan data (exon microarray summarized to genes)

# read in expression data (rows = genes, columns = samples)
expDf <- read.csv("../data/BrainSpan/array_expression_matrix.processed.csv",
	header=T,stringsAsFactor=F)

# subset to data from samples in the four frontal cortex regions
expDf <- expDf %>% select(contains("gene"),
	contains("MFC"),contains("DFC"),contains("OFC"),contains("VFC"))

# read in gene and sample meta data
geneDf <- read.csv("../data/BrainSpan/rows_metadata.filtered.csv", stringsAsFactor=F)
sampDf <- read.csv("../data/BrainSpan/columns_metadata.csv", stringsAsFactors=F)
# only keep samples from cortex regions
sampDf <- subset(sampDf,Identifier %in% colnames(expDf)) 

# read in age -> developmental stage/time point mapping
# then assign each sample to a time point
ageDf <- read.csv("../data/BrainSpan/age.mapping.csv", stringsAsFactors=F)
sampDf$TimePoint <- plyr::mapvalues(sampDf$age, ageDf$Age, ageDf$Developmental.Period)

# store number of samples for each time point
sampCounts <- as.data.frame(table(sampDf$TimePoint))
colnames(sampCounts) <- c('TimePoint','NumSamples')


# ----------------------------------------------------------------------------------------
# read in gene sets

# ASD genes (Satterstrom)
asdDf <- read.table('../data/Satterstrom2020_TableS2_Autosomal.txt',header=T,sep='\t')
# columns to use: hugoGene, qval_dnccPTV


# COMBINED PPI network
intDf <- read.table('../data/ASD_MasterInteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)

baits <- sort(unique(intDf$Bait))
baits <- baits[baits != 'COMBINED']

intDf <- subset(intDf, Dataset=='COMBINED')[,c('Dataset','Gene','IsInteractor')]


# Satterstrom 25% FDR genes excluding baits
s25_noBaits <- asdDf$hugoGene[asdDf$qval_dnccPTV <= 0.25 & !asdDf$hugoGene %in% baits]

# Satterstrom 25% FDR genes in COMBINED network, exclude baits
s25_netGenes <- intersect(s25_noBaits,subset(intDf,Dataset=='COMBINED' & IsInteractor)$Gene)


# named list of lists containing each gene set
geneSets <- list(
	'5% FWER' = asdDf$hugoGene[asdDf$qval_dnccPTV <= 1.4e-4],
	'25% FDR' = asdDf$hugoGene[asdDf$qval_dnccPTV <= 0.25],
	
	'Network' = subset(intDf,Dataset=='COMBINED' & IsInteractor)$Gene,

	'25% FDR Network' = s25_netGenes,
	
	'Random' = c())
	

# filter gene sets containing only genes with available BrainSpan data
filteredSets <- lapply(geneSets, function(x) {x[x %in% geneDf$gene_symbol]})


# random genes (matching # of genes in '25% FDR Network')
set.seed(1234)
filteredSets$'Random' <- sample(geneDf$gene_symbol,length(filteredSets$'25% FDR Network'))

# check # of genes in each gene set pre- and post-filter
geneCounts <- merge(plyr::ldply(geneSets, .fun=length, .id='Set'),
	plyr::ldply(filteredSets, .fun=length, .id='Set'),by='Set')
colnames(geneCounts) <- c('Set','PreGenes','NumGenes')
geneCounts <- geneCounts[order(geneCounts$Set),]
rownames(geneCounts) <- NULL

geneCounts


# ----------------------------------------------------------------------------------------
# extract expression values for each gene set at each time point

# iterate through each gene set
meltDf <- NULL
for (s in names(filteredSets)) {
	tempDf <- expDf[geneDf$row_num[geneDf$gene_symbol %in% filteredSets[[s]]],]
	tempDf <- melt(tempDf,id.vars='gene')
	tempDf$Set <- s
	
	meltDf <- rbind(meltDf,tempDf)
}

meltDf$TimePoint <- plyr::mapvalues(meltDf$variable, sampDf$Identifier, sampDf$TimePoint)


# ----------------------------------------------------------------------------------------
# calculate median, SD, and SE expression of each gene set at each time point

# median
tpMedian <- melt(dcast(meltDf, Set~TimePoint, value.var='value', fun.aggregate=median), 
	id.vars='Set', variable.name='TimePoint', value.name='Median')

# standard deviation
tpSd <- melt(dcast(meltDf, Set~TimePoint, value.var='value', fun.aggregate=sd),
	id.vars='Set', variable.name='TimePoint', value.name='SD')

# number of data points (# genes x # samples) used at each time point
tpCounts <- melt(dcast(meltDf, Set~TimePoint, value.var='value', fun.aggregate=length),
	id.vars='Set', variable.name='TimePoint', value.name='N')

# merge stats into one data frame
mergeDf <- merge(merge(tpMedian, tpSd, by=c('Set','TimePoint')), tpCounts, by=c('Set','TimePoint'))

# calculate SE: standard deviation divided by sqrt(number of data points)
mergeDf$SE <- mergeDf$SD/sqrt(mergeDf$N)

# add number of genes and samples as columns
mergeDf$nSamples <- plyr::mapvalues(mergeDf$TimePoint,sampCounts$TimePoint,sampCounts$NumSamples)
mergeDf$nGenes <- plyr::mapvalues(mergeDf$Set,geneCounts$Set,geneCounts$NumGenes)


# order table by time points (chronologically) and gene sets
mergeDf$Set <- factor(mergeDf$Set, levels=names(geneSets))
mergeDf$TimePoint <- factor(mergeDf$TimePoint, levels=unique(ageDf$Developmental.Period))
mergeDf <- mergeDf[order(mergeDf$Set,mergeDf$TimePoint),]


# ----------------------------------------------------------------------------------------
# generate BrainSpan plot for Satterstrom 25% FDR vs. Network genes vs. their overlap

# data to for line plot
lineSets <- c('5% FWER','25% FDR','Network','25% FDR Network','Random')
lineDf <- subset(mergeDf,Set %in% lineSets)
lineDf$Set <- factor(lineDf$Set, levels=lineSets)


pdf("../output/ASD_SatterstromOverlap_BrainSpanPlot.pdf",height=4,width=5)

ggplot(lineDf,aes(x=TimePoint,y=Median)) +

# plot median +- SE for gene sets
geom_line(aes(color=Set,group=Set),size=0.5) + geom_point(aes(color=Set,group=Set)) +
geom_errorbar(aes(ymin=Median-SE,ymax=Median+SE,color=Set,group=Set),size=0.5,width=0.2) +
scale_color_manual(name='Gene set',
	labels=c('5% FWER','25% FDR','Network','25% FDR + Network Overlap','Random'),
	values=c('black','darkorange2','purple','magenta','aquamarine3')) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=11,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (\u00B1 SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(color=guide_legend(nrow=2,byrow=T))

dev.off()



##########################################################################################
# generate BrainSpan plots for:
# (1) COMBINED network genes
# (2) iN vs brain interactors of SCN2A, SYNGAP1, SHANK3

# read in expresion values (computed by Kalliopi Tsafou using same methods as above)
df <- read.table('../data/BrainSpan/ASD_iN_brain_BrainSpanValues.tsv',
	header=T,sep='\t',stringsAsFactors=F)

# columns:
# variable=developmental stage
# gene=gene set
# value.x=median expression value of all genes in each developmental stage in each set
# value.y=standard deviation for all genes set in each developmental stage in each set
# value=total number of genes expressed in all patients per developmental stage per set
# Freq=number of samples per developmental stage

df$variable <- factor(df$variable,
	levels=c('Early prenatal','Early mid-prenatal','Late mid-prenatal','Late prenatal',
		'Early infancy','Late infancy','Early childhood','Late childhood',
		'Adolescence','Adulthood'))


# ----------------------------------------------------------------------------------------
# COMBINED network BrainSpan plot

# data to for line plot
geneSets <- c('Baits','iNs','FWER <= 0.05','Random set')
lineDf <- subset(df,gene %in% geneSets)
lineDf$gene <- factor(lineDf$gene,levels=geneSets)

# data for shading (Satterstrom gene sets)
fdrSets <- c('FWER <= 0.05','FDR <= 0.1','FDR <= 0.25')
fdrDf <- subset(df,gene %in% fdrSets)
fdrDf$gene <- factor(fdrDf$gene,levels=fdrSets)

fdrDf$fdr <- 1 # FDR grouping for color scale
fdrDf$fdr[fdrDf$gene=='FDR <= 0.1'] <- 2
fdrDf$fdr[fdrDf$gene=='FDR <= 0.25'] <- 3

fdrDf$ribbon.min <- 6
fdrDf$ribbon.min[fdrDf$gene=='FWER <= 0.05'] <- fdrDf$value.x[fdrDf$gene=='FDR <= 0.1']
fdrDf$ribbon.min[fdrDf$gene=='FDR <= 0.1'] <- fdrDf$value.x[fdrDf$gene=='FDR <= 0.25']
# code above requires table is sorted so time points for each gene set are in same order


pdf('../output/ASD_BrainSpanPlots.pdf',height=4,width=5)

ggplot(lineDf,aes(x=variable,y=value.x)) + 

# plot Satterstrom gene sets as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=variable,ymax=value.x,ymin=ribbon.min,
	fill=gene,group=gene),linetype='blank') +
scale_fill_manual(name='Satterstrom',
	labels=c('FWER <= 0.05','FDR <= 0.1','FDR <= 0.25'),
	values=c('grey60','grey80','grey100')) +

# plot median +- SE for index, network, FWER <= 0.05 Satterstrom, random genes
geom_line(aes(color=gene,group=gene),size=0.5) + geom_point(aes(color=gene,group=gene)) +
geom_errorbar(aes(ymin=value.x-(value.y/sqrt(value)),ymax=value.x+(value.y/sqrt(value)),
	color=gene,group=gene),size=0.5,width=0.2) +
scale_color_manual(name='Gene set',
	labels=c('Index','Network','Satterstrom','Random'),
	values=c('red2','purple','black','aquamarine3')) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=11,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=1,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=2))
	

# ----------------------------------------------------------------------------------------
# iN vs brain interactors BrainSpan plot

# data to for line plot
geneSets <- c('iNs_unique','Cortex_unique','iNs_cortex_common','FWER <= 0.05','Random set')
lineDf <- subset(df,gene %in% geneSets)
lineDf$gene <- factor(lineDf$gene,levels=geneSets)

# data for shading (Satterstrom gene sets)
fdrSets <- c('FWER <= 0.05','FDR <= 0.1','FDR <= 0.25')
fdrDf <- subset(df,gene %in% fdrSets)
fdrDf$gene <- factor(fdrDf$gene,levels=fdrSets)

fdrDf$fdr <- 1 # FDR grouping for color scale
fdrDf$fdr[fdrDf$gene=='FDR <= 0.1'] <- 2
fdrDf$fdr[fdrDf$gene=='FDR <= 0.25'] <- 3

fdrDf$ribbon.min <- 6
fdrDf$ribbon.min[fdrDf$gene=='FWER <= 0.05'] <- fdrDf$value.x[fdrDf$gene=='FDR <= 0.1']
fdrDf$ribbon.min[fdrDf$gene=='FDR <= 0.1'] <- fdrDf$value.x[fdrDf$gene=='FDR <= 0.25']
# code above requires table is sorted so time points for each gene set are in same order


# PLOT
ggplot(lineDf,aes(x=variable,y=value.x)) + 

# plot Satterstrom gene sets as shaded regions
geom_ribbon(data=fdrDf,mapping=aes(x=variable,ymax=value.x,ymin=ribbon.min,
	fill=gene,group=gene),linetype='blank') +
scale_fill_manual(name='Satterstrom',
	labels=c('FWER <= 0.05','FDR <= 0.1','FDR <= 0.25'),
	values=c('grey60','grey80','grey100')) +

# plot median +- SE for gene sets
geom_line(aes(color=gene,group=gene),size=0.5) + geom_point(aes(color=gene,group=gene)) +
geom_errorbar(aes(ymin=value.x-(value.y/sqrt(value)),ymax=value.x+(value.y/sqrt(value)),
	color=gene,group=gene),size=0.5,width=0.2) +
scale_color_manual(name='Gene set',
	labels=c('iN','Brain','Shared','Satterstrom','Random'),
	values=c('purple','darkorange2','blue','black','aquamarine3')) +

# pre vs. postnatal line
geom_vline(xintercept=4.5,linetype='dashed') +
annotate("text",x=4.6,y=11,hjust=0,size=4,label='Birth') +

xlab("Developmental stage") +
ylab("Median expression (± SE)") +
theme_classic() +
theme(legend.position='bottom',legend.box='vertical',legend.margin=margin(),
	legend.title=element_text(size=8),legend.text=element_text(size=8),
	axis.text.x=element_text(angle=45,hjust=1)) +
guides(fill=guide_legend(order=1,override.aes=list(linetype=1,color='black')),
	color=guide_legend(order=2))
	
dev.off()

