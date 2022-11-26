##########################################################################################
## GTEx tissue enrichment analysis
##
## Author: Kalliopi Tsafou, Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------------------------------
# networks to test: COMBINED vs. InWeb-based network

intDf <- read.table('../data/ASD_MasterInteractorTable_withInWeb.txt',
	header=T,sep='\t',stringsAsFactors=F)

# list of baits (for removing them from background population in GTEx below)
baits <- unique(intDf$Bait[intDf$Bait!='COMBINED'])
length(baits) # 13

intDf <- subset(intDf,(Dataset=='COMBINED' | Dataset=='InWeb') & IsInteractor)
table(intDf$Dataset) # check number of genes in each network

# ----------------------------------------------------------------------------------------
# read in GTEx gene tissue specificity scores and define top 10% of genes in each tissue
# to be tissue-specific genes
# GTEx data from: https://alkesgroup.broadinstitute.org/LDSCORE/LDSC_SEG_ldscores/tstats/

# ENSG ID to gene name mapping
ensgDf <- read.table('../data/Homo_sapiens.GRCh38.84.GeneAnnotations.txt',header=T,sep='\t')

# data frame to store all GTEx sets
otherDf <- NULL

# GTEx data for all tissues
gtDf <- read.table('../data/GTEx/GTEx.tstat.tsv',header=T,sep='\t',row.names=1)
gtNameDf <- read.table('../data/GTEx/GTEx.tissue.names.tsv',header=T,sep='\t')

gtDf$gene <- ensgDf$GeneName[match(rownames(gtDf),ensgDf$GeneID)]
gtDf <- subset(gtDf,!is.na(gene) & !duplicated(gene) & !gene %in% baits)
dim(gtDf) # 24083 unique gene names x 54 columns (gene name + 53 tissues)

top10_cutoff <- floor(nrow(gtDf)*.1) # number of genes in top 10% (2408)

for (tissue in names(gtDf)[names(gtDf)!='gene']) {
	# save top 10% of genes ranked by dereacing t score
	tissueGenes <- gtDf$gene[order(gtDf[,tissue],decreasing=T)][1:top10_cutoff]
	tissueName <- gtNameDf$TissueName[gtNameDf$Tissue==tissue]
	otherDf <- rbind(otherDf,data.frame(SetType='GTEx_all',SetName=tissueName,gene=tissueGenes))
}

# GTEx data for brain tissues
gtbDf <- read.table('../data/GTEx/GTEx_brain.tstat.tsv',header=T,sep='\t',row.names=1)
gtbNameDf <- read.table('../data/GTEx/GTEx_brain.tissue.names.tsv',header=T,sep='\t')

gtbDf$gene <- ensgDf$GeneName[match(rownames(gtbDf),ensgDf$GeneID)]
gtbDf <- subset(gtbDf,!is.na(gene) & !duplicated(gene) & !gene %in% baits)
dim(gtbDf) # 24083 unique gene names x 14 columns (gene name + 13 brain tissues)
# same # of genes as for all tissues above

for (tissue in names(gtbDf)[names(gtbDf)!='gene']) {
	# save top 10% of genes ranked by dereacing t score
	tissueGenes <- gtbDf$gene[order(gtbDf[,tissue],decreasing=T)][1:top10_cutoff]
	tissueName <- gtbNameDf$TissueName[gtbNameDf$Tissue==tissue]
	otherDf <- rbind(otherDf,data.frame(SetType='GTEx_brain',SetName=tissueName,gene=tissueGenes))
}


# ----------------------------------------------------------------------------------------
# perform overlap enrichment tests

# results table
outDf <- NULL

for (setType in unique(otherDf$SetType)) {
	setDf <- subset(otherDf,SetType==setType)

	for (setName in unique(setDf$SetName)) {
	
		for (network in unique(intDf$Dataset)) {
			
			popGenes <- gtDf$gene

			# success: network genes in population
			intGenes <- subset(intDf,Dataset==network & Gene %in% popGenes)$Gene

			# sample: tissue-specific genes in population
			setGenes <- subset(setDf,SetName==setName & gene %in% popGenes)$gene
			
			# success in sample: overlap between network and tissue-specific genes
			overlapGenes <- sort(intersect(intGenes,setGenes))

			m <- length(intGenes)
			n <- length(setdiff(popGenes,intGenes))
			k <- length(setGenes)
			q <- length(overlapGenes)

			p <- phyper(q-1,m,n,k,lower.tail=F)

			outDf <- rbind(outDf,data.frame(SetType=setType,SetName=setName,Network=network,
				PopCount=length(popGenes),IntCount=m,SetCount=k,OverlapCount=q,OverlapP=p,
				OverlapGenes=paste(overlapGenes,collapse=',')))

		}
	}
}

# output results table
write.table(outDf,'../output/ASD_GtexEnrichment.txt',sep='\t',row.names=F,quote=F)


# ----------------------------------------------------------------------------------------
# plot all tissue enrichment

pdf ('../output/ASD_GtexEnrichment_BarPlots.pdf',height=3,width=3)

for (net in c('COMBINED','InWeb')) {
	df <- subset(outDf,SetType=='GTEx_all' & Network==net)
	df$Category <- gtNameDf$TissueCategory[match(df$SetName,gtNameDf$TissueName)]
	nSets <- nrow(df)

	# set tissue plotting order, by tissue category, then by p-value
	df <- df[order(df$Category,df$OverlapP),]
	df$SetName <- factor(df$SetName,levels=rev(df$SetName))
	
	p <- ggplot(df,aes(x=-log10(OverlapP),y=SetName,fill=Category)) + geom_col() +

		geom_vline(xintercept=-log10(0.05),linetype='dashed',size=0.25) +
		geom_vline(xintercept=-log10(0.05/nSets),linetype='dashed',size=0.25) +

		xlab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
		ylab('Tissue') +
		labs(fill='Tissue category') +
		ggtitle(paste(net,'network')) +
		
		theme_bw() +
		theme(axis.title=element_text(size=8),
			axis.text.y=element_blank(),axis.text.x=element_text(size=7),
			legend.title=element_text(size=8),legend.text=element_text(size=7),
			legend.key.size=unit(0.25,'lines'))

	print(p)

}


# ----------------------------------------------------------------------------------------
# plot brain tissue enrichment

for (net in c('COMBINED','InWeb')) {
	df <- subset(outDf,SetType=='GTEx_brain' & Network==net)
	nSets <- nrow(df)
	df$Sig <- ifelse(df$OverlapP < 0.05,'Nominal','None')
	df$Sig[df$OverlapP < 0.05/nSets] <- 'Bonferroni'

	df$SigText <- ifelse(df$OverlapP < 0.05,df$OverlapCount,'')

	# set tissue plotting order, by p-value
	df <- df[order(df$OverlapP),]
	df$SetName <- factor(df$SetName,levels=rev(df$SetName))

	p <- ggplot(df,aes(x=-log10(OverlapP),y=SetName,fill=Sig)) + geom_col() +

		geom_vline(xintercept=-log10(0.05),linetype='dashed',size=0.25) +
		geom_vline(xintercept=-log10(0.05/nSets),linetype='dashed',size=0.25) +

		geom_text(aes(label=SigText),hjust=-0.1,size=3) +

		scale_fill_manual(values = c('None'='black','Nominal'='sienna1','Bonferroni'='red')) +
		scale_x_continuous(expand=expansion(mult=c(0,0.2),add=c(0,0))) +

		xlab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
		ylab('Brain tissue') +
		ggtitle(paste(net,'network')) +
		theme_bw() + theme(legend.position='none',
			axis.title=element_text(size=8),axis.text=element_text(size=7))

	print(p)
}

dev.off()

