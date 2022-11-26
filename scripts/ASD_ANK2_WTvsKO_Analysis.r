##########################################################################################
## compare IP-MS results from ANK2-WT vs ANK2-KO (giant ANK2 KO; or ANK2-MUT)
## neural progenitor cells
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)
library(ggrepel)
library(genoppi)

set.seed(1234)

# ----------------------------------------------------------------------------------------
# read in Genoppi analysis results and define significant interactors

# WT vs IgG results
wtDf <- read.table('../data/ANK2_WTvsMUT/ANK2-WT.pep1.bait2.MinImp.GenoppiStats.txt',
	header=T,sep='\t',stringsAsFactors=F)
wtDf$significant <- wtDf$logFC > 0 & wtDf$FDR <= 0.1

# KO vs IgG results
mutDf <- read.table('../data/ANK2_WTvsMUT/ANK2-MUT.pep1.bait2.MinImp.GenoppiStats.txt',
	header=T,sep='\t',stringsAsFactors=F)
mutDf$significant <- mutDf$logFC > 0 & mutDf$FDR <= 0.1

# WT vs KO results
vsDf <- read.table('../data/ANK2_WTvsMUT/ANK2-WTvsMUT.pep1.bait2.MinImp.GenoppiStats.txt',
	header=T,sep='\t',stringsAsFactors=F)

# annotate proteins based on significance in WT vs IgG or KO vs IgG results
vsDf$sig_wt <- vsDf$gene %in% wtDf$gene[wtDf$significant]
vsDf$sig_mut <- vsDf$gene %in% mutDf$gene[mutDf$significant]

# define WT or KO-specific interactors
vsDf$sig_group <- 'None'
vsDf$sig_group[vsDf$sig_wt & !vsDf$sig_mut] <- 'WT'
vsDf$sig_group[vsDf$sig_mut & !vsDf$sig_wt] <- 'KO'
vsDf$sig_group <- factor(vsDf$sig_group,levels=c('None','WT','KO'))

wtDf$sig_group <- 'None'
wtDf$sig_group[wtDf$gene %in% vsDf$gene[vsDf$sig_group=='WT']] <- 'WT'
wtDf$sig_group <- factor(wtDf$sig_group,levels=c('None','WT'))

mutDf$sig_group <- 'None'
mutDf$sig_group[mutDf$gene %in% vsDf$gene[vsDf$sig_group=='KO']] <- 'KO'
mutDf$sig_group <- factor(mutDf$sig_group,levels=c('None','KO'))

# proteins analyzed in both WT vs IgG and KO vs IgG experiments
combDf <- merge(wtDf,mutDf,by='accession_number')
combDf$sig_group <- as.character(combDf$sig_group.x)
combDf$sig_group[combDf$sig_group.y=='KO'] <- 'KO'
combDf$sig_group <- factor(combDf$sig_group,levels=c('None','WT','KO'))

# alternative combDf that assign artificial values to proteins only identified in one experiment
combDf_NA <- merge(wtDf,mutDf,by='accession_number',all=T)
combDf_NA$sig_group <- vsDf$sig_group[match(combDf_NA$accession_number,vsDf$accession_number)]
combDf_NA$gene <- vsDf$gene[match(combDf_NA$accession_number,vsDf$accession_number)]

#minVal <- min(combDf_NA$logFC.x,combDf_NA$logFC.y,na.rm=T) # -5.387135
combDf_NA$logFC.x[is.na(combDf_NA$logFC.x)] <- -6
combDf_NA$logFC.y[is.na(combDf_NA$logFC.y)] <- -6


# ----------------------------------------------------------------------------------------
# volcano plots with specificity annotations

bait <- 'ANK2'

pdf('../output/ASD_ANK2_WTvsKO_Plots.pdf',height=5.5,width=5)

# WT vs KO volcano: color based on WT or KO-specificity
ggplot(vsDf,aes(x=logFC,y=-log10(pvalue),color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote("WT/KO "*log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='Specificity',values=c('grey','red','blue')) + 

# label bait
geom_point(subset(vsDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color='orange') +
geom_point(subset(vsDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),size=2,color="black",shape=1) +

# label proteins that are WT or KO specific
geom_point(subset(vsDf, sig_group=='WT' | sig_group=='KO'),
	mapping=aes(x=logFC, y=-log10(pvalue),color=sig_group),size=2) +
geom_point(subset(vsDf, sig_group=='WT' | sig_group=='KO'),
	mapping=aes(x=logFC, y=-log10(pvalue)),size=2,color="black",shape=1) +

geom_text_repel(subset(vsDf, gene==bait | sig_group=='WT' | sig_group=='KO'),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2, max.overlaps=1000) +

ggtitle('ANK2 IP in WT vs. KO cells') +
theme_bw() + theme(axis.line=element_line(color="grey"),legend.position='bottom',
	legend.box="vertical", legend.margin=margin())


# WT vs IgG volcano: color based on WT-specificty
ggplot(wtDf,aes(x=logFC,y=-log10(pvalue),color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote("WT/IgG "*log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='WT-Specific',labels=c('No','Yes'),values=c('grey','red')) + 

# label bait
geom_point(subset(wtDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color='orange') +
geom_point(subset(wtDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color="black",shape=1) +

# label proteins that are WT specific
geom_point(subset(wtDf, sig_group=='WT'),mapping=aes(x=logFC, y=-log10(pvalue),
	color=sig_group),size=2) +
geom_point(subset(wtDf, sig_group=='WT'),mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color="black",shape=1) +

geom_text_repel(subset(wtDf, gene==bait | sig_group=='WT'),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2, max.overlaps=1000) +

ggtitle('ANK2 vs IgG IP in WT cells') +
theme_bw() + theme(axis.line=element_line(color="grey"),legend.position='bottom',
	legend.box="vertical", legend.margin=margin())


# KO vs IgG volcano: color based on KO-specificty
ggplot(mutDf,aes(x=logFC,y=-log10(pvalue),color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote("KO/IgG "*log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='KO-Specific',labels=c('No','Yes'),values=c('grey','blue')) + 

# label bait
geom_point(subset(mutDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color='orange') +
geom_point(subset(mutDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color="black",shape=1) +

# label proteins that are KO specific
geom_point(subset(mutDf, sig_group=='KO'),mapping=aes(x=logFC, y=-log10(pvalue),
	color=sig_group),size=2) +
geom_point(subset(mutDf, sig_group=='KO'),mapping=aes(x=logFC, y=-log10(pvalue)),size=2,color="black",shape=1) +

geom_text_repel(subset(mutDf, gene==bait | sig_group=='KO'),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2, max.overlaps=1000) +

ggtitle('ANK2 vs IgG IP in KO cells') +
theme_bw() + theme(axis.line=element_line(color="grey"),legend.position='bottom',
	legend.box="vertical", legend.margin=margin())


# WT/IgG vs KO/IgG log2 FC scatter plot: color based on WT or KO-specificity
ggplot(combDf,aes(x=logFC.x,y=logFC.y,color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote("WT/IgG "*log[2]*"(Fold change)")) +
ylab(bquote("KO/IgG "*log[2]*"(Fold change)")) +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='Specificity',values=c('grey','red','blue')) + 

# label bait
geom_point(subset(combDf, gene.x==bait), mapping=aes(x=logFC.x, y=logFC.y),
	size=2,color='orange') +
geom_point(subset(combDf, gene.x==bait), mapping=aes(x=logFC.x, y=logFC.y),
	size=2,color="black",shape=1) +

# label proteins that are WT or KO specific
geom_point(subset(combDf, sig_group=='WT' | sig_group=='KO'),
	mapping=aes(x=logFC.x, y=logFC.y,color=sig_group),size=2) +
geom_point(subset(combDf, sig_group=='WT' | sig_group=='KO'),
	mapping=aes(x=logFC.x, y=logFC.y),size=2,color="black",shape=1) +

geom_text_repel(subset(combDf, gene.x==bait | sig_group=='WT' | sig_group=='KO'),
	mapping=aes(label=gene.x),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2, max.overlaps=1000) +

ggtitle('ANK2 WT/IgG vs KO/IgG log2 FC') +
theme_bw() + theme(axis.line=element_line(color="grey"),legend.position='bottom',
	legend.box="vertical", legend.margin=margin())


# WT/IgG vs KO/IgG log2 FC scatter plot: color based on WT or KO-specificity
# plot proteins that were analyzed in only one experiment using artificial min values
ggplot(combDf_NA,aes(x=logFC.x,y=logFC.y,color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote("WT/IgG "*log[2]*"(Fold change)")) +
ylab(bquote("KO/IgG "*log[2]*"(Fold change)")) +

# dashed lines to highlight artificial min logFC
geom_hline(yintercept=-5.5, color="black",linetype='dashed') +
geom_vline(xintercept=-5.5, color="black",linetype='dashed') +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='Specificity',values=c('grey','red','blue')) + 

# label bait
geom_point(subset(combDf_NA, gene.x==bait), mapping=aes(x=logFC.x, y=logFC.y),
	size=2,color='orange') +
geom_point(subset(combDf_NA, gene.x==bait), mapping=aes(x=logFC.x, y=logFC.y),
	size=2,color="black",shape=1) +

# label proteins that are WT or KO specific
geom_point(subset(combDf_NA, sig_group=='WT' | sig_group=='KO'),
	mapping=aes(x=logFC.x, y=logFC.y,color=sig_group),size=2) +
geom_point(subset(combDf_NA, sig_group=='WT' | sig_group=='KO'),
	mapping=aes(x=logFC.x, y=logFC.y),size=2,color="black",shape=1) +

geom_text_repel(subset(combDf_NA, gene==bait | sig_group=='WT' | sig_group=='KO'),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2, max.overlaps=1000) +

ggtitle('ANK2 WT/IgG vs KO/IgG log2 FC') +
theme_bw() + theme(axis.line=element_line(color="grey"),legend.position='bottom',
	legend.box="vertical", legend.margin=margin())

dev.off()


# ----------------------------------------------------------------------------------------
# GO cellular component enrichment analysis

overlapDf <- NULL
overlapDf_genome <- NULL

ccList <- unique(goa_cc_table$GO.name) # use GO terms stored in genoppi R package

for (cc in ccList) {

	ccDf <- data.frame(gene=subset(goa_cc_table,GO.name==cc)$Gene.symbol,significant=TRUE)

	# hypergeometric test for WT-specific ints
	intDf <- data.frame(gene=vsDf$gene,significant=vsDf$sig_group=='WT')

	results <- calc_hyper(intDf,ccDf,data.frame(intersectN=F),bait='ANK2')
	
	resDf <- data.frame(Specificity='WT',GO_CC=cc,
		results$statistics[,c('population_count','success_count','sample_count',
			'successInSample_count','pvalue')])
			
	colnames(resDf) <- c('Specificity','GO_CC','Population',
		'SpecificInts','CC','Overlap','Pvalue')
	resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
	overlapDf <- rbind(overlapDf,resDf)

	# hypergeometric test for KO-specific ints
	intDf <- data.frame(gene=vsDf$gene,significant=vsDf$sig_group=='KO')

	results <- calc_hyper(intDf,ccDf,data.frame(intersectN=F),bait='ANK2')
	
	resDf <- data.frame(Specificity='KO',GO_CC=cc,
		results$statistics[,c('population_count','success_count','sample_count',
			'successInSample_count','pvalue')])
			
	colnames(resDf) <- c('Specificity','GO_CC','Population',
		'SpecificInts','CC','Overlap','Pvalue')
	resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
	overlapDf <- rbind(overlapDf,resDf)

}

# output enrichment results
write.table(overlapDf,'../output/ASD_ANK2_WTvsKO_CellularComponentEnrichment.txt',
	sep='\t',quote=F,row.names=F)

