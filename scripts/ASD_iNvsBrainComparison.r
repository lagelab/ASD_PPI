##########################################################################################
## Compare SCN2A, SYANK3, SYNGAP1 interactors found in induced neurons (iN) vs brain IPs
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(readxl)
library(ggplot2)
library(ggrepel)
library(genoppi)


# ----------------------------------------------------------------------------------------
# read in iN interactors and non-interactors

baits <- c('SCN2A','SHANK3','SYNGAP1')

inDf <- read.table('../data/ASD_MasterInteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)
inDf <- subset(inDf, Dataset %in% baits)
inDf$gene <- inDf$Gene
inDf$significant <- inDf$IsInteractor


# add number of IPs for each interactor
numIpDf <- read.table('../data/ASD_IntSubsetTable.txt',
	header=T,sep='\t',stringsAsFactors=F)
inDf$NumIPs <- NA
for (i in 1:nrow(inDf)) {
	if (inDf$significant[i]) {
		inDf$NumIPs[i] <- nrow(subset(numIpDf,Bait==inDf$Bait[i] & 
			Interactor==inDf$gene[i]))
	}
}


# combined iN ints and non-ints of the 3 baits
inInts <- setdiff(unique(subset(inDf,significant)$gene),baits)
inNonInts <- setdiff(setdiff(unique(subset(inDf,!significant)$gene),inInts),baits)


# ----------------------------------------------------------------------------------------
# read in brain interactors, assess overlap with iN interactors, plot overlaid volcanos

overlapDf <- NULL
brainInts <- NULL
brainNonInts <- NULL

pdf('../output/ASD_BrainVolcanos_iNoverlay.pdf',height=4,width=4)

for (bait in baits) {

	# Brain ints and non-ints
	brainDf <- read_excel("../data/ASD_GenoppiResults.xlsx",
		sheet=paste(bait,'_Brain',sep=''))
		
	brainDf$logFC <- as.numeric(brainDf$logFC)
	brainDf$pvalue <- as.numeric(brainDf$pvalue)
	brainDf$FDR <- as.numeric(brainDf$FDR)
	
	brainDf$significant <- brainDf$logFC > 0 & brainDf$FDR <= 0.1
	

	# overlap with iN ints
	brainDf$InNeuron <- brainDf$gene %in% subset(inDf,Bait==bait & significant)$gene


	# (1) volcano plot overlaid with iN interactors
	p <- ggplot(brainDf,aes(x=logFC,y=-log10(pvalue))) +
	
		# zero lines on both axes
		geom_hline(yintercept=0,color="grey") + geom_vline(xintercept=0,color="grey") +
		
		# axis labels
		xlab(bquote(log[2]*"(Fold change)")) +
		ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
		
		# plot all proteins (green = significant, blue = not significant)
		geom_point(alpha=0.5, size=1.5,
			color=ifelse(brainDf$significant,"springgreen3","royalblue2")) +
		
		# label bait
		geom_point(subset(brainDf,gene==bait & significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +
		
		# label iN interactors (yellow = significant, white = not sig)
		geom_point(subset(brainDf,InNeuron & significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="yellow") +
		geom_point(subset(brainDf,InNeuron & !significant),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="white") +
			
		geom_point(subset(brainDf,gene==bait | InNeuron),
			mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="black",shape=1) +
			
		# bait name label
		geom_text_repel(subset(brainDf,gene==bait),
			mapping=aes(label=gene),arrow=arrow(length=unit(0.015,'npc')),
			box.padding=unit(0.15,"lines"),point.padding=unit(0.2,"lines"),
			color="black",size=3,max.overlaps=1000) +

		ggtitle(paste(bait,'_Brain',sep='')) + theme_bw() + 
		theme(axis.line=element_line(color="grey"))

	print(p)
	
	
	# (2) iN vs brain interactor overlap enrichment test
	brainDf <- brainDf[,c('gene','significant')]

	# add filtered non-interactors
	# (i.e. detected proteins filtered out prior to Genoppi analysis)
	filterNonInts <- read.table(paste('../data/BrainIPs_Filtered/',bait,
		'.pep1.bait2.FilteredNonInts.txt',sep=''),
		header=T,sep='\t',stringsAsFactors=F)$gene
	
	brainDf <- rbind(brainDf,data.frame(gene=filterNonInts,significant=FALSE))

	# hypergeometric test to compare iN vs. brain ints (intersectN=T)
	results <- calc_hyper(subset(inDf,Bait==bait)[,c('gene','significant')],
		data.frame(listName=bait,brainDf),
		data.frame(listName=bait,intersectN=T),bait=bait)
		
	resDf <- data.frame(Bait=bait,
		NeuronInt=length(unique(subset(inDf,Bait==bait & significant)$gene)),
		BrainInt=length(unique(subset(brainDf,significant & gene!=bait)$gene)),
		results$statistics[,c('population_count','success_count',
			'sample_count','successInSample_count','pvalue')])
			
	colnames(resDf) <- c('Bait','NeuronInt','BrainInt',
		'Population','NeuronInPop','BrainInPop','Overlap','Pvalue')
	resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
	overlapDf <- rbind(overlapDf,resDf)


	# hypergeometric test for iN ints detected in only 1 IP
	results <- calc_hyper(subset(inDf,Bait==bait &
		(NumIPs==1 | is.na(NumIPs)))[,c('gene','significant')],
		data.frame(listName=bait,brainDf),
		data.frame(listName=bait,intersectN=T),bait=bait)
		
	resDf <- data.frame(Bait=paste(bait,', =1IP',sep=''),
		NeuronInt=length(unique(subset(inDf,Bait==bait & significant & NumIPs==1)$gene)),
		BrainInt=length(unique(subset(brainDf,significant & gene!=bait)$gene)),
		results$statistics[,c('population_count','success_count',
			'sample_count','successInSample_count','pvalue')])
			
	colnames(resDf) <- c('Bait','NeuronInt','BrainInt',
		'Population','NeuronInPop','BrainInPop','Overlap','Pvalue')
	resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
	overlapDf <- rbind(overlapDf,resDf)

	
	# hypergeometric test for iN ints detected in > 1 IP
	results <- calc_hyper(subset(inDf,Bait==bait &
		(NumIPs>1 | is.na(NumIPs)))[,c('gene','significant')],
		data.frame(listName=bait,brainDf),
		data.frame(listName=bait,intersectN=T),bait=bait)
		
	resDf <- data.frame(Bait=paste(bait,', >1IP',sep=''),
		NeuronInt=length(unique(subset(inDf,Dataset==bait & significant & NumIPs>1)$gene)),
		BrainInt=length(unique(subset(brainDf,significant & gene!=bait)$gene)),
		results$statistics[,c('population_count','success_count',
			'sample_count','successInSample_count','pvalue')])
	
	colnames(resDf) <- c('Bait','NeuronInt','BrainInt',
		'Population','NeuronInPop','BrainInPop','Overlap','Pvalue')
	resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
	overlapDf <- rbind(overlapDf,resDf)
	

	# add brain ints and non-ints to combined lists
	brainInts <- c(brainInts,subset(brainDf,significant)$gene)
	brainNonInts <- c(brainNonInts,subset(brainDf,!significant)$gene)
}

dev.off()


# combined brain ints and non-ints of the 3 baits
brainInts <- setdiff(unique(brainInts),baits)
brainNonInts <- setdiff(setdiff(unique(brainNonInts),brainInts),baits)

inAllDf <- data.frame(gene=c(inInts,inNonInts),
	significant=c(rep(TRUE,length(inInts)),rep(FALSE,length(inNonInts))))
brainAllDf <- data.frame(gene=c(brainInts,brainNonInts),
	significant=c(rep(TRUE,length(brainInts)),rep(FALSE,length(brainNonInts))))


# hypergeometric test for combined ints
results <- calc_hyper(inAllDf,data.frame(listName='COMBINED',brainAllDf),
	data.frame(listName='COMBINED',intersectN=T))
resDf <- data.frame(Bait=paste(baits,collapse='_'),
	NeuronInt=length(inInts),BrainInt=length(brainInts),
	results$statistics[,c('population_count','success_count',
		'sample_count','successInSample_count','pvalue')])
colnames(resDf) <- c('Bait','NeuronInt','BrainInt',
	'Population','NeuronInPop','BrainInPop','Overlap','Pvalue')
resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
overlapDf <- rbind(overlapDf,resDf)


# hypergeometric test for iN ints detected in only 1 IP
inInts <- setdiff(unique(subset(inDf,significant & NumIPs==1)$gene),baits)
inNonInts <- setdiff(setdiff(unique(subset(inDf,!significant)$gene),inInts),baits)
inAllDf <- data.frame(gene=c(inInts,inNonInts),
	significant=c(rep(TRUE,length(inInts)),rep(FALSE,length(inNonInts))))
results <- calc_hyper(inAllDf,data.frame(listName='COMBINED',brainAllDf),
	data.frame(listName='COMBINED',intersectN=T))
resDf <- data.frame(Bait=paste(paste(baits,collapse='_'),', =1IP',sep=''),
	NeuronInt=length(inInts),BrainInt=length(brainInts),
	results$statistics[,c('population_count','success_count',
		'sample_count','successInSample_count','pvalue')])
colnames(resDf) <- c('Bait','NeuronInt','BrainInt',
	'Population','NeuronInPop','BrainInPop','Overlap','Pvalue')
resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
overlapDf <- rbind(overlapDf,resDf)


# hypergeometric test for iN ints detected in > 1 IP
inInts <- setdiff(unique(subset(inDf,significant & NumIPs>1)$gene),baits)
inNonInts <- setdiff(setdiff(unique(subset(inDf,!significant)$gene),inInts),baits)
inAllDf <- data.frame(gene=c(inInts,inNonInts),
	significant=c(rep(TRUE,length(inInts)),rep(FALSE,length(inNonInts))))
results <- calc_hyper(inAllDf,data.frame(listName='COMBINED',brainAllDf),
	data.frame(listName='COMBINED',intersectN=T))
resDf <- data.frame(Bait=paste(paste(baits,collapse='_'),', >1IP',sep=''),
	NeuronInt=length(inInts),BrainInt=length(brainInts),
	results$statistics[,c('population_count','success_count',
		'sample_count','successInSample_count','pvalue')])
colnames(resDf) <- c('Bait','NeuronInt','BrainInt',
	'Population','NeuronInPop','BrainInPop','Overlap','Pvalue')
resDf$OverlapGenes <- paste(results$genes[[1]]$successInSample_genes,collapse=',')
overlapDf <- rbind(overlapDf,resDf)


# output overlap enrichment test results
write.table(overlapDf,'../output/ASD_iNvsBrain_OverlapEnrichment.txt',
	sep='\t',quote=F,row.names=F)

