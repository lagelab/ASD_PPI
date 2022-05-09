##########################################################################################
## Generate volcano and scatter plots for ASD IP-MS datasets
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(readxl)
library(ggplot2)
library(ggrepel)
library(genoppi)

# ----------------------------------------------------------------------------------------
# read in Genoppi analysis results each IP

df <- NULL
sheetList <- excel_sheets("../data/ASD_GenoppiResults.xlsx")

for (i in 2:27) { # only want the 26 IPs done in neurons
	tempDf <- read_excel("../data/ASD_GenoppiResults.xlsx",sheet=sheetList[i])
	tempDf$IP <- sheetList[i]

	# index protein	
	bait <- strsplit(sheetList[i],'_')[[1]][1]
	tempDf$bait <- bait
	
	# define significant proteins
	tempDf$significant = tempDf$logFC>0 & tempDf$FDR<=0.1
	
	# InWeb interactors
	tempInwebDf <- get_inweb_list(bait)
	if (!is.null(tempInwebDf)) {
		tempDf$IsInWeb <- tempDf$gene %in% subset(tempInwebDf,significant)$gene
	} else {
		tempDf$IsInWeb <- FALSE
	}
	
	# save data to master tables for plotting
	df <- rbind(df,tempDf)
}


# ----------------------------------------------------------------------------------------
# Scatter and InWeb volcano plots for 13 representative IPs

ipList <- c('ADNP_CNCR1','ANK2_WH','ARID1B_CNCR1','CHD8_CNCR2','CTNNB1_CNCR1',
	'DYRK1A_Broad','GIGYF1_CNCR2','MED13L_CNCR1','PTEN_WH1','SCN2A_Broad',
	'SHANK3_Broad','SYNGAP1_WH1','TLK2_CNCR1')

pdf("../output/ASD_Scatter_InWebVolcanos_13IPs.pdf",height=4,width=4)

for (ip in ipList) {

	subDf <- subset(df,IP==ip)

	r <- cor(subDf$rep1,subDf$rep2)

	# SCATTER
	p1 <- ggplot(subDf,aes(x=rep1,y=rep2)) +

		# plot all proteins (green = significant, blue = not significant)
		geom_point(alpha=0.5,size=1.5,
			color=ifelse(subDf$significant,"springgreen3","royalblue2")) +

		# label bait (red = significant, orange = not significant)
		geom_point(subset(subDf,gene==bait & significant),
			mapping=aes(x=rep1,y=rep2),size=2,color="red") +
		geom_point(subset(subDf,gene==bait & !significant),
			mapping=aes(x=rep1,y=rep2),size=2,color="orange") +
			
		geom_point(subset(subDf,gene==bait),mapping=aes(x=rep1,y=rep2),
			size=2,color="black",shape=1) +
		
		geom_text_repel(subset(subDf,gene==bait),mapping=aes(label=gene),
			arrow=arrow(length=unit(0.015,'npc')),box.padding=unit(0.15,"lines"),
			point.padding=unit(0.2,"lines"),color="black",size=3) +

		# identity line, title (with correlation), theme
		geom_abline(intercept=0, slope=1, linetype="longdash", size=0.2) +
		ggtitle(paste(ip,", correlation:",format(r,digits=3),sep='')) +
		theme_bw() + theme(axis.line=element_line(color="grey"))

	print(p1)


	# INWEB VOLCANO
	p2 <- ggplot(subDf,aes(x=logFC,y=-log10(pvalue))) +
	
		# zero lines on both axes
		geom_hline(yintercept=0,color="grey") + geom_vline(xintercept=0,color="grey") +
		
		# axis labels
		xlab(bquote(log[2]*"(Fold change)")) +
		ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
		
		# plot all proteins (green = significant, blue = not significant)
		geom_point(alpha=0.5, size=1.5,
			color=ifelse(subDf$significant,"springgreen3","royalblue2")) +
		
		# label bait (red = significant, orange = not significant)
		geom_point(subset(subDf,gene==bait & significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +
		geom_point(subset(subDf,gene==bait & !significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="orange") +
		
		# label InWeb interactors (yellow = significant, white = not sig)
		geom_point(subset(subDf,IsInWeb & significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="yellow") +
		geom_point(subset(subDf,IsInWeb & !significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="white") +
			
		geom_point(subset(subDf,gene==bait | IsInWeb),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="black",shape=1) +
			
		geom_text_repel(subset(subDf,gene==bait | IsInWeb),
			mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
			box.padding=unit(0.15,"lines"),point.padding=unit(0.2,"lines"),
			color="black",size=3,max.overlaps=1000) +

		ggtitle(ip) + theme_bw() + 
		theme(axis.line=element_line(color="grey"))

	print(p2)
	
}

dev.off()


# ----------------------------------------------------------------------------------------
# ANK2 annotated volcano plot

pdf("../output/ASD_ANK2_PTEN_DYRK1A_Volcanos.pdf",height=4,width=4)


subDf <- subset(df,IP=='ANK2_WH')
orangeGenes <- c('ACTB','GAP43','C3','CACNA2D1')


ggplot(subDf, aes(x=logFC, y=-log10(pvalue))) +

# zero lines on both axes
geom_hline(yintercept=0,color="grey") + geom_vline(xintercept=0,color="grey") +

# axis labels
xlab(bquote(log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

# plot all proteins (green = significant, blue = not significant)
geom_point(alpha=0.5, size=1.5,
	color=ifelse(subDf$significant,"springgreen3","royalblue2")) +
	
# label bait
geom_point(subset(subDf,gene==bait & significant),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +

# label InWeb interactors (yellow = significant, white = not sig)
		geom_point(subset(subDf,IsInWeb & significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="yellow") +
		geom_point(subset(subDf,IsInWeb & !significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="white") +
			
# highlight genes
geom_point(subset(subDf, gene %in% orangeGenes),
	mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +

geom_point(subset(subDf, gene==bait | IsInWeb | gene %in% orangeGenes),
	mapping=aes(x=logFC, y=-log10(pvalue)),size=2, color="black", shape=1) +

geom_text_repel(subset(subDf, gene==bait | IsInWeb | gene %in% orangeGenes),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
	box.padding=unit(0.15,"lines"),point.padding=unit(0.2,"lines"),color="black",size=3) +

# title and theme
ggtitle('ANK2_WH IP-MS data') +
theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))


# ----------------------------------------------------------------------------------------
# PTEN annotated volcano plot

subDf <- subset(df,IP=='PTEN_WH1')
orangeGenes <- c('AKAP8L')


ggplot(subDf, aes(x=logFC, y=-log10(pvalue))) +

# zero lines on both axes
geom_hline(yintercept=0,color="grey") + geom_vline(xintercept=0,color="grey") +

# axis labels
xlab(bquote(log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

# plot all proteins (green = significant, blue = not significant)
geom_point(alpha=0.5, size=1.5,
	color=ifelse(subDf$significant,"springgreen3","royalblue2")) +
	
# label bait
geom_point(subset(subDf,gene==bait & significant),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +

# highlight genes
geom_point(subset(subDf, gene %in% orangeGenes),
	mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +

geom_point(subset(subDf, gene==bait | gene %in% orangeGenes),
	mapping=aes(x=logFC, y=-log10(pvalue)),size=2, color="black", shape=1) +

geom_text_repel(subset(subDf, gene==bait | gene %in% orangeGenes),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
	box.padding=unit(0.15,"lines"),point.padding=unit(0.2,"lines"),color="black",size=3) +

# title and theme
ggtitle('PTEN_WH1 IP-MS data') +
theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))


# ----------------------------------------------------------------------------------------
# DYRK1A annotated volcano plot

subDf <- subset(df,IP=='DYRK1A_Broad')
orangeGenes <- c('IGF2BP1','IGF2BP2','IGF2BP3')


ggplot(subDf, aes(x=logFC, y=-log10(pvalue))) +

# zero lines on both axes
geom_hline(yintercept=0,color="grey") + geom_vline(xintercept=0,color="grey") +

# axis labels
xlab(bquote(log[2]*"(Fold change)")) +
ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +

# plot all proteins (green = significant, blue = not significant)
geom_point(alpha=0.5, size=1.5,
	color=ifelse(subDf$significant,"springgreen3","royalblue2")) +
	
# label bait
geom_point(subset(subDf,gene==bait & significant),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +

# highlight genes
geom_point(subset(subDf, gene %in% orangeGenes),
	mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +

geom_point(subset(subDf, gene==bait | gene %in% orangeGenes),
	mapping=aes(x=logFC, y=-log10(pvalue)),size=2, color="black", shape=1) +

geom_text_repel(subset(subDf, gene==bait | gene %in% orangeGenes),
	mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
	box.padding=unit(0.15,"lines"),point.padding=unit(0.2,"lines"),color="black",size=3) +

# title and theme
ggtitle('DYRK1A_Broad IP-MS data') +
theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))

dev.off()