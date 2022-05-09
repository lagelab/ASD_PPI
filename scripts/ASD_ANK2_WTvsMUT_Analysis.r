##########################################################################################
## compare IP-MS results from ANK2-WT vs ANK2-MUT (giant ANK2 KO) neural progenitor cells
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)
library(ggrepel)
library(genoppi)


# ----------------------------------------------------------------------------------------
# read in Genoppi analysis results and define significant interactors

# WT vs IgG results
wtDf <- read.table('../data/ANK2_WTvsMUT/ANK2-WT.pep1.bait2.MinImp.GenoppiStats.txt',
	header=T,sep='\t',stringsAsFactors=F)
wtDf$significant <- wtDf$logFC > 0 & wtDf$FDR <= 0.1

# MUT vs IgG results
mutDf <- read.table('../data/ANK2_WTvsMUT/ANK2-MUT.pep1.bait2.MinImp.GenoppiStats.txt',
	header=T,sep='\t',stringsAsFactors=F)
mutDf$significant <- mutDf$logFC > 0 & mutDf$FDR <= 0.1

# WT vs MUT results
vsDf <- read.table('../data/ANK2_WTvsMUT/ANK2-WTvsMUT.pep1.bait2.MinImp.GenoppiStats.txt',
	header=T,sep='\t',stringsAsFactors=F)

# annotate proteins based on significance in WT or MUT vs IgG results
vsDf$sig_wt <- vsDf$gene %in% wtDf$gene[wtDf$significant]
vsDf$sig_mut <- vsDf$gene %in% mutDf$gene[mutDf$significant]

# define WT or MUT-specific interactors
vsDf$sig_group <- 'None'
vsDf$sig_group[vsDf$sig_wt & !vsDf$sig_mut] <- 'WT'
vsDf$sig_group[vsDf$sig_mut & !vsDf$sig_wt] <- 'MUT'
vsDf$sig_group <- factor(vsDf$sig_group,levels=c('None','WT','MUT'))

wtDf$sig_group <- 'None'
wtDf$sig_group[wtDf$gene %in% vsDf$gene[vsDf$sig_group=='WT']] <- 'WT'
wtDf$sig_group <- factor(wtDf$sig_group,levels=c('None','WT'))

mutDf$sig_group <- 'None'
mutDf$sig_group[mutDf$gene %in% vsDf$gene[vsDf$sig_group=='MUT']] <- 'MUT'
mutDf$sig_group <- factor(mutDf$sig_group,levels=c('None','MUT'))


# ----------------------------------------------------------------------------------------
# volcano plots with specificity annotations

bait <- 'ANK2'
plotDf <- vsDf

pdf('../output/ASD_ANK2_WTvsMUT_Volcanos.pdf',height=5.5,width=5)

# WT vs MUT volcano: color based on WT or MUT-specificity
ggplot(plotDf,aes(x=logFC,y=-log10(pvalue),color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote(log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='Specificity',values=c('grey','red','blue')) + 

# label bait
geom_point(subset(plotDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color='orange') +
geom_point(subset(plotDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),size=2,color="black",shape=1) +

# label proteins that are WT or MUT specific
geom_point(subset(plotDf, sig_group=='WT' | sig_group=='MUT'),
	mapping=aes(x=logFC, y=-log10(pvalue),color=sig_group),size=2) +
geom_point(subset(plotDf, sig_group=='WT' | sig_group=='MUT'),
	mapping=aes(x=logFC, y=-log10(pvalue)),size=2,color="black",shape=1) +

geom_text_repel(subset(plotDf, gene==bait | sig_group=='WT' | sig_group=='MUT'),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2) +

ggtitle('ANK2 IP in WT vs. MUT cells') +
theme_bw() + theme(axis.line=element_line(color="grey"),legend.position='bottom',
	legend.box="vertical", legend.margin=margin())


# WT vs IgG volcano: color based on WT-specificty
ggplot(wtDf,aes(x=logFC,y=-log10(pvalue),color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote(log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='WT-Specific',labels=c('No','Yes'),values=c('grey','red')) + 

# label bait
geom_point(subset(wtDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color='orange') +
geom_point(subset(wtDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color="black",shape=1) +

# label proteins that are WT or MUT specific
geom_point(subset(wtDf, sig_group=='WT'),mapping=aes(x=logFC, y=-log10(pvalue),
	color=sig_group),size=2) +
geom_point(subset(wtDf, sig_group=='WT'),mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color="black",shape=1) +

geom_text_repel(subset(wtDf, gene==bait | sig_group=='WT'),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2) +

ggtitle('ANK2 vs IgG IP in WT cells') +
theme_bw() + theme(axis.line=element_line(color="grey"),legend.position='bottom',
	legend.box="vertical", legend.margin=margin())


# MUT vs IgG volcano: color based on MUT-specificty
ggplot(mutDf,aes(x=logFC,y=-log10(pvalue),color=sig_group)) +
geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="grey") +
xlab(bquote(log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

geom_point(size=1.5,alpha=0.7) +
scale_color_manual(name='MUT-Specific',labels=c('No','Yes'),values=c('grey','blue')) + 

# label bait
geom_point(subset(mutDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color='orange') +
geom_point(subset(mutDf, gene==bait), mapping=aes(x=logFC, y=-log10(pvalue)),
	size=2,color="black",shape=1) +

# label proteins that are WT or MUT specific
geom_point(subset(mutDf, sig_group=='MUT'),mapping=aes(x=logFC, y=-log10(pvalue),
	color=sig_group),size=2) +
geom_point(subset(mutDf, sig_group=='MUT'),mapping=aes(x=logFC, y=-log10(pvalue)),size=2,color="black",shape=1) +

geom_text_repel(subset(mutDf, gene==bait | sig_group=='MUT'),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015, 'npc')),
	box.padding=unit(0.15, "lines"), point.padding=unit(0.2, "lines"),
	color="black", size=2) +

ggtitle('ANK2 vs IgG IP in MUT cells') +
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

	# hypergeometric test for MUT-specific ints
	intDf <- data.frame(gene=vsDf$gene,significant=vsDf$sig_group=='MUT')

	results <- calc_hyper(intDf,ccDf,data.frame(intersectN=F),bait='ANK2')
	
	resDf <- data.frame(Specificity='MUT',GO_CC=cc,
		results$statistics[,c('population_count','success_count','sample_count',
			'successInSample_count','pvalue')])
			
	colnames(resDf) <- c('Specificity','GO_CC','Population',
		'SpecificInts','CC','Overlap','Pvalue')
	resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
	overlapDf <- rbind(overlapDf,resDf)

}

# output enrichment results
write.table(overlapDf,'../output/ASD_ANK2_WTvsMUT_CellularComponentEnrichment.txt',
	sep='\t',quote=F,row.names=F)

