##########################################################################################
## Compare and plot co-expression z-scores between bait-interactor pairs vs. other
## gene pairs, using co-expression data derived from 4 transcriptomic datasets:
## Stickels2021, Maynard2021, Velmeshev2019, BrainSpan
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(genoppi)
library(ggpubr)

# ----------------------------------------------------------------------------------------
# read in bait-gene pairs in which genes are (1) interactors, (2) non-interactors,
# or (3) InWeb interactors

tempDf <- read.table('../data/ASD_MasterInteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)
tempDf <- subset(tempDf,!grepl('COMBINED',Dataset) & !grepl('_',Dataset))

# (1) bait-interactor pairs
intDf <- data.frame(subset(tempDf,IsInteractor)[c('Dataset','Bait','Gene')],Type='Int')

# (2) bait-non-interactor pairs
nonDf <- data.frame(subset(tempDf,!IsInteractor)[c('Dataset','Bait','Gene')],Type='NonInt')


baits <- unique(intDf$Bait)

# (3) bait-InWeb interactor pairs
webDf <- NULL
for (b in baits) {
	x <- subset(get_inweb_list(b),significant)$gene
	webDf <- rbind(webDf,
		data.frame(Dataset=b,Bait=b,Gene=x,Type='InWeb'))
}

# combine dfs and add COMBINED Dataset
mergeDf <- rbind(intDf,nonDf,webDf)

combDf <- mergeDf
combDf$Dataset <- 'COMBINED'
mergeDf <- rbind(mergeDf,combDf)


# ----------------------------------------------------------------------------------------
# read in co-expression data

coexpFiles <- c('../data/CoExp/Stickels2021','../data/CoExp/Maynard2021',
	'../data/CoExp/Velmeshev2019','../data/CoExp/BrainSpan')
metrics <- c('FisherNegLogP','FisherNegLogP','PropRho','PropRho')
titles <- c('Stickels 2021','Maynard 2021','Velmeshev 2019','BrainSpan')

plotDf <- NULL

for (i in 1:length(titles)) {
	print(titles[i])
	
	inFile <- gzfile(paste(coexpFiles[i],'_ASD-Baits_CoExp',metrics[i],'.txt.gz',sep=''),'rt')
	
	metricDf <- read.table(inFile,header=T,sep='\t',
		stringsAsFactors=F,row.names=1,check.names=F)
		
	close(inFile)

	if (metrics[i]=='FisherNegLogP') {
		# fix any infinite numbers outputted by Python
		rowNames <- rownames(metricDf)
		metricDf[metricDf=='inf'] <- Inf
		metricDf <- data.frame(apply(metricDf,2,function(x) as.numeric(x)))
		rownames(metricDf) <- rowNames

		# set Inf values to just above max of non-Inf value (per bait/row)
		metricDf <- data.frame(t(apply(metricDf,1,function(x) {
			ifelse(is.infinite(x),ceiling(max(x[!is.infinite(x)])),x)})))
	}

	# rank-based inverse normal transformation to get z-scores (per bait/row)
	metricDf <- data.frame(t(apply(metricDf,1,function(x) {
		qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))})))

	# save z-scores for bait-gene pairs involving Int, NonInt, InWeb genes
	iDf <- subset(mergeDf, Bait %in% rownames(metricDf) & Gene %in% names(metricDf))
	matchRows <- match(iDf$Bait,rownames(metricDf))
	matchCols <- match(iDf$Gene,names(metricDf))
	for (k in 1:length(matchRows)) {
		iDf$Metric[k] <- metricDf[matchRows[k],matchCols[k]]
	}

	# save z-scores for bait-gene pairs for All genes
	tempDf <- NULL
	for (b in unique(iDf$Bait)) {
		tempDf <- rbind(tempDf,data.frame(Dataset=b,Bait=b,
			Gene=names(metricDf),Type='All',Metric=as.numeric(metricDf[b,])))
	}
	iDf <- rbind(iDf,tempDf)

	# save z-scores for COMBINED Dataset for All genes
	tempDf$Dataset <- 'COMBINED'
	iDf <- rbind(iDf,tempDf)

	# filter out bait-bait pairs
	iDf <- subset(iDf,Bait!=Gene)
	plotDf <- rbind(plotDf,data.frame(iDf,DataSource=titles[i]))

}


# ----------------------------------------------------------------------------------------
# plot and compare co-expression z-scores of Int, NonInt, InWeb, and All genes

pdf('../output/ASD_CoExp_ViolinPlots.pdf',height=2.5,width=7.5)

# set plotting order
plotDf$Type <- factor(plotDf$Type,levels=unique(plotDf$Type))

plotDf$DataSource <- factor(plotDf$DataSource,levels=titles)

dOrder <- sort(unique(plotDf$Dataset))
plotDf$Dataset <- factor(plotDf$Dataset,levels=c('COMBINED',dOrder[dOrder!='COMBINED']))

for (d in levels(plotDf$Dataset)) {

	# set Wilcoxon test parameters
	varPairs <- combn(levels(plotDf$Type),2)
	testList <- NULL
	for (col in 1:ncol(varPairs)) { testList[[col]] <- varPairs[,col] }
	pCutoffs <- c(0,0.05/length(testList),0.05,1)
	pSymbols <- c('**','*','ns')

	# plot score distributions
	plot <- ggviolin(subset(plotDf,Dataset==d),x='Type',y='Metric',color='Type') +
		geom_boxplot(aes(color=Type),width=0.2,outlier.size=0.2) +
		
		facet_wrap(~DataSource,ncol=4,scales='free_y') +

		# pairwise Wilcoxon tests
		stat_compare_means(comparisons=testList,method='wilcox.test',
			symnum.args=list(cutpoints=pCutoffs,symbols=pSymbols),
			hide.ns=T,size=3,vjust=0.5) +

		# label size of each group
		stat_summary(aes(label=..y..,color=Type,
			y=stage(Metric,after_stat=min(plotDf$Metric))),
			fun=length,geom='text',size=2.5,vjust=1.4,
			position=position_dodge(width=0.75)) +

		scale_color_brewer(palette='Set1') +
		scale_y_continuous(expand=expansion(mult=c(0.1,0.05),add=c(0,0))) +
		ggtitle(d) +
		xlab('Gene type') + ylab(bquote("Co-expression"~italic(Z)*"-score")) +
		theme_bw() + theme(legend.position='none',axis.text=element_text(size=8))
	
	print(plot)

}

dev.off()


# ----------------------------------------------------------------------------------------
# representative plot: Maynard2021 COMBINED results

# set Wilcoxon test parameters
testList <- list(c('Int','NonInt'),c('Int','InWeb'),c('Int','All'))
# correct for all 6 pairs of comparison even though only showing 3 in plot
pCutoffs <- c(0,0.05/6,0.05,1)
pSymbols <- c('**','*','ns')


pdf('../output/ASD_CoExp-Maynard2021_ViolinPlot.pdf',height=2,width=2.5)

# plot score distributions
ggviolin(subset(plotDf,Dataset=='COMBINED' & DataSource=='Maynard 2021'),
	x='Type',y='Metric',color='Type')+
geom_boxplot(aes(color=Type),width=0.2,outlier.size=0.2) +

# pairwise Wilcoxon tests
stat_compare_means(comparisons=testList,method='wilcox.test',
	symnum.args=list(cutpoints=pCutoffs,symbols=pSymbols),hide.ns=T,size=3,vjust=0.5) +

# label size of each group
stat_summary(aes(label=..y..,color=Type,y=stage(Metric,after_stat=min(plotDf$Metric))),
fun=length,geom='text',size=2.5,vjust=0,position=position_dodge(width=0.75)) +

scale_color_brewer(palette='Set1') +
scale_y_continuous(expand=expansion(mult=c(0.1,0.05),add=c(0,0))) +
xlab('Gene type') + ylab(bquote("Co-expression"~italic(Z)*"-score")) +
theme_bw() + theme(legend.position='none',axis.text=element_text(size=8),
	axis.title=element_text(size=9))

dev.off()

