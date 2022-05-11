##########################################################################################
## Enrichment analysis for IGF2BP1-3 targets
## (1) overlap enrichment between IGF2BP1-3 targets vs. Satterstrom ASD genes
## (2) plot ASD common variant enrichment of IGF2BP1-3 targets (from MAGMA)
## (3) overlap enrichment between IGF2BP1-3 targets vs. ASD PPI networks
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)


# ----------------------------------------------------------------------------------------
# read in Satterstrom ASD genes
asdDf <- read.table('../data/Satterstrom2020_TableS2_Autosomal.txt',
	header=T,sep='\t',stringsAsFactors=F)


# read in ints/non-ints in ASD PPI networks
intDf <- read.table('../data/ASD_MasterInteractorTable.txt',header=T,stringsAsFactors=F)
intDf <- subset(intDf,!grepl('_',Dataset))

# number of interactor genes per network
netSizeDf <- data.frame(table(subset(intDf, IsInteractor)$Dataset))
colnames(netSizeDf) <- c('Dataset','IntCount')
	

# read in IGF2BP1-3 RNA targets
targDf <- read.table('../data/Huang2018_IGF2BP1-3_RNAtargets.txt',
	header=T,stringsAsFactors=F) # Protein_Gene. RNA_Gene

# number of target genes per IGF2BP (or union of IGF2BP1-3)
sizeDf <- data.frame(table(targDf$Protein_Gene))
colnames(sizeDf) <- c('IGF2BP','TargetCount')
sizeDf <- rbind(sizeDf,
	data.frame(IGF2BP='IGF2BP1-3',TargetCount=length(unique(targDf$RNA_Gene))))


# ----------------------------------------------------------------------------------------
# overlap enrichment between Satterstrom ASD genes vs. IGF2BP targets

asdOverlapDf <- NULL

pop <- asdDf$gene # autosomal genes with q-values in Satterstrom paper
sample <- subset(asdDf,qval_dnccPTV<=0.1)$gene
notSample <- subset(asdDf,qval_dnccPTV>0.1)$gene

# ASD genes vs. targets of each IGF2BP
for (i in unique(targDf$Protein_Gene)) {

	success <- subset(targDf,Protein_Gene==i)$RNA_Gene
	success <- success[success %in% pop]
	successInSample <- sort(sample[sample %in% success])
	hyperP <- phyper(length(successInSample)-1,length(sample),
		length(notSample),length(success),lower.tail=F)

	overlapBaits <- successInSample[successInSample %in% unique(intDf$Bait)]

	asdOverlapDf <- rbind(asdOverlapDf,data.frame(List1="ASD_Satterstrom",
		List2=paste(i),OverlapCount=length(successInSample),
		List1Count=length(sample),List2Count=length(success),
		PopCount=length(pop),HyperP=hyperP,
		OverlapGenes=paste(successInSample,collapse=','),
		OverlapBaits=paste(overlapBaits,collapse=',')))
}

# ASD genes vs. union of IGF2BP targets
success <- unique(targDf$RNA_Gene)
success <- success[success %in% pop]
successInSample <- sort(sample[sample %in% success])
hyperP <- phyper(length(successInSample)-1,length(sample),
	length(notSample),length(success),lower.tail=F)

overlapBaits <- successInSample[successInSample %in% unique(intDf$Bait)]

asdOverlapDf <- rbind(asdOverlapDf,data.frame(List1="ASD_Satterstrom",
	List2=paste("IGF2BP1-3"),OverlapCount=length(successInSample),
	List1Count=length(sample),List2Count=length(success),
	PopCount=length(pop),HyperP=hyperP,
	OverlapGenes=paste(successInSample,collapse=','),
	OverlapBaits=paste(overlapBaits,collapse=',')))

# output results table
write.table(asdOverlapDf,'../output/ASD_IGF2BPtargets-ASDgenes_OverlapEnrichment.txt',
	row.names=F,quote=F,sep='\t')


# ----------------------------------------------------------------------------------------
# plot IGF2BP targets rare (above) and common (MAGMA) variant enrichment results

# RARE VARIANTS
df <- data.frame(IGF2BP=asdOverlapDf$List2,
	Overlap=asdOverlapDf$OverlapCount,P=asdOverlapDf$HyperP)

# set target list label and plotting order
df$IGF2BP <- factor(df$IGF2BP,levels=c('IGF2BP3','IGF2BP2','IGF2BP1','IGF2BP1-3'))
df <- df[order(df$IGF2BP),]
df$IGF2BP <- paste(df$IGF2BP,'\n(',
	sizeDf$TargetCount[match(df$IGF2BP,sizeDf$IGF2BP)],')',sep='')
df$IGF2BP <- factor(df$IGF2BP,levels=df$IGF2BP)

# set p-value text color/font based on nominal or Bonferroni significance
# Bonferroni: adjust for 4 IGF2BP target lists
df$Sig <- ifelse(df$P < 0.05,"Nominal","None")
df$Sig[df$P<0.05/4] <- "Bonferroni"
df$Sig <- factor(df$Sig, levels=c("None","Nominal","Bonferroni"))


# bar plot
pdf("../output/ASD_IGF2BPtargets-ASDgenes_OverlapEnrichmentPlot.pdf",height=2,width=2)

ggplot(df,aes(x=IGF2BP,y=-log10(P),fill=Sig)) + geom_col(show.legend=F) +

geom_text(aes(label=Overlap),hjust=-0.1,size=3) + 

geom_hline(yintercept=-log10(0.05/4), linetype="dashed",size=0.25) +
scale_fill_manual(values=c("None"="black","Nominal"="sienna1","Bonferroni"="red")) +
scale_y_continuous(expand=expansion(mult=c(0,0.2),add=c(0,0))) +

xlab('RNA targets') + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
ggtitle('rare variant') +
coord_flip() + theme_bw() +
theme(plot.title=element_text(size=9),
	axis.title=element_text(size=9),axis.text=element_text(size=8))


# COMMON VARIANTS
df <- read.table('../data/ASD_IGF2BPtargets_MagmaResults.txt',
	header=T,sep='\t',stringsAsFactors=F)

df <- subset(df,TRAIT=='ASD')

# set target list label and plotting order
df$IGF2BP <- df$VARIABLE
df$IGF2BP[df$IGF2BP=='COMBINED'] <- 'IGF2BP1-3'
df$IGF2BP <- factor(df$IGF2BP,levels=c('IGF2BP3','IGF2BP2','IGF2BP1','IGF2BP1-3'))
df <- df[order(df$IGF2BP),]
df$IGF2BP <- paste(df$IGF2BP,'\n(',
	sizeDf$TargetCount[match(df$IGF2BP,sizeDf$IGF2BP)],')',sep='')
df$IGF2BP <- factor(df$IGF2BP,levels=df$IGF2BP)

# set p-value text color/font based on nominal or Bonferroni significance
# Bonferroni: adjust for 4 IGF2BP target lists
df$Sig <- ifelse(df$P < 0.05,"Nominal","None")
df$Sig[df$P<0.05/4] <- "Bonferroni"
df$Sig <- factor(df$Sig, levels=c("None","Nominal","Bonferroni"))


# bar plot
ggplot(df,aes(x=IGF2BP,y=-log10(P),fill=Sig)) + geom_col(show.legend=F) +

geom_text(aes(label=formatC(BETA,digits=2)),hjust=-0.1,size=3) + 

geom_hline(yintercept=-log10(0.05/4), linetype="dashed",size=0.25) +

scale_fill_manual(values=c("None"="black","Nominal"="sienna1","Bonferroni"="red")) +
scale_y_continuous(expand=expansion(mult=c(0,0.5),add=c(0,0))) +

xlab('RNA targets') + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)"))  + ggtitle('common variant') +
coord_flip() + theme_bw() + 
theme(plot.title=element_text(size=9),
	axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()


# ----------------------------------------------------------------------------------------
# overlap enrichment between ASD PPI network interactors vs. IGF2BP targets

overlapDf <- NULL
for (d in unique(intDf$Dataset)) {

	tempDf <- subset(intDf,Dataset==d)
	pop <- tempDf$Gene
	sample <- subset(tempDf,IsInteractor)$Gene
	notSample <- subset(tempDf,!IsInteractor)$Gene
	
	# each IGF2BP targets
	for (i in unique(targDf$Protein_Gene)) {
		
		success <- subset(targDf,Protein_Gene==i)$RNA_Gene
		success <- success[success %in% pop]
		successInSample <- sort(sample[sample %in% success])
		hyperP <- phyper(length(successInSample)-1,length(sample),
			length(notSample),length(success),lower.tail=F)

		overlapDf <- rbind(overlapDf,data.frame(Dataset=d,Protein=i,
			OverlapCount=length(successInSample),TargetCount=length(success),
			IntCount=length(sample),PopCount=length(pop),HyperP=hyperP,
			OverlapGenes=paste(successInSample,collapse=',')))
	}
	
	# union of IGF2BP1-3 targets
	success <- unique(targDf$RNA_Gene)
	success <- success[success %in% pop]
	successInSample <- sort(sample[sample %in% success])
	hyperP <- phyper(length(successInSample)-1,length(sample),
		length(notSample),length(success),lower.tail=F)

	overlapDf <- rbind(overlapDf,data.frame(Dataset=d,Protein="IGF2BP1-3",
		OverlapCount=length(successInSample),TargetCount=length(success),
		IntCount=length(sample),PopCount=length(pop),HyperP=hyperP,
		OverlapGenes=paste(successInSample,collapse=',')))
	
}

# output results table
write.table(overlapDf,'../output/ASD_IGF2BPtargets-Networks_OverlapEnrichment.txt',
	row.names=F,quote=F,sep='\t')


# ----------------------------------------------------------------------------------------
# plot IGF2BP targets vs. PPI networks enrichment results

df <- data.frame(IGF2BP=overlapDf$Protein,Network=overlapDf$Dataset,
	Overlap=overlapDf$OverlapCount,P=overlapDf$HyperP)
	
# set target list label and plotting order
df$IGF2BP <- factor(df$IGF2BP,levels=c('IGF2BP3','IGF2BP2','IGF2BP1','IGF2BP1-3'))
df <- df[order(df$IGF2BP),]
df$IGF2BP <- paste(df$IGF2BP,'\n(',
	sizeDf$TargetCount[match(df$IGF2BP,sizeDf$IGF2BP)],')',sep='')
df$IGF2BP <- factor(df$IGF2BP,levels=unique(df$IGF2BP))

# set network label and plotting order
baitOrder <- sort(unique(overlapDf$Dataset))
baitOrder <- c('COMBINED',baitOrder[baitOrder!='COMBINED'])
df$Network <- factor(df$Network,levels=baitOrder)
df <- df[order(df$Network),]
df$Network <- paste(df$Network,'\n(',
	netSizeDf$IntCount[match(df$Network,netSizeDf$Dataset)],')',sep='')
df$Network <- factor(df$Network,levels=unique(df$Network))
 
# set significance marker text for nominal for Bonferroni significance
# Bonferroni: adjust for 14 networks and 4 IGF2BP target lists
df$Sig <- ifelse(df$P < 0.05,"*","")
df$Sig[df$P < 0.05/4/14] <- "**"

df$SigText <- ifelse(df$P < 0.05,paste(df$Overlap,df$Sig,sep=''),'')

# show baits that interact with any IGF2BPs in blue
goodNets <- c('COMBINED','CHD8','DYRK1A','SCN2A','SHANK3','SYNGAP1','ANK2','TLK2')
baitColors <- ifelse(rowSums(sapply(goodNets,grepl,levels(df$Network)))==1,
	'blue','grey30')

# heat map
pdf("../output/ASD_IGF2BPtargets-Networks_OverlapEnrichmentPlot.pdf",height=2.5,width=5.5)

ggplot(df,aes(x=Network,y=IGF2BP,fill=-log10(P))) +
geom_tile() + geom_text(aes(label=SigText),size=3) +

scale_fill_gradient(name=bquote(-log[10]*"("*italic(.("P"))*"-value)"),
	low="white",high="red") +

xlab('Network') + ylab('RNA targets') +
theme_bw() +
theme(legend.position="bottom",legend.key.size=unit(0.7,"line"),
	legend.title=element_text(size=9),legend.text=element_text(size=8),
	legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,10,-5,-10),
	axis.title=element_text(size=9),axis.text.y=element_text(size=8),
	axis.text.x=element_text(color=baitColors,size=8,angle=90,vjust=0.5,hjust=1))

dev.off()

