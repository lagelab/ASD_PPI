##########################################################################################
## Rare variant risk enrichment analysis of ASD PPI networks
## using gene-based scores and KS test
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(genoppi) # for pLI scores in gnomad_table
library(ggplot2)


# ----------------------------------------------------------------------------------------
# gene-based scores to use in KS tests

# gnomAD pLI scores
pliDf <- unique(subset(gnomad_table,!is.na(gnomad_table$pLI)))
scoreDf <- data.frame(variable='gnomAD pLI',gene=pliDf$gene,score=pliDf$pLI)

# ASD q-values (Satterstrom et al.; columns to use: hugoGene, qval_dnccPTV)
asdDf <- read.table('../data/Satterstrom2020_TableS2_Autosomal.txt',header=T,sep='\t')
scoreDf <- rbind(scoreDf,data.frame(variable='ASD Q-value',
	gene=asdDf$hugoGene,score=asdDf$qval_dnccPTV))

# DD p-values (Kaplanis et al.; columns to use: symbol, denovoWEST_p_full)
ddDf <- read.table('../data/Kaplanis2019_TableS2.txt',header=T,sep='\t')
ddDf <- subset(ddDf,!is.na(denovoWEST_p_full))
scoreDf <- rbind(scoreDf,data.frame(variable='DD P-value',
	gene=ddDf$symbol,score=ddDf$denovoWEST_p_full))
        
# SCZ/SCHEMA p-values (Singh et al.; columns to use: Gene.Symbol, P.meta, Q.meta)
schemaDf <- read.csv('../data/SCHEMA2020_TableS5.csv')
schemaDf <- subset(schemaDf,!is.na(P.meta))
scoreDf <- rbind(scoreDf,data.frame(variable='SCHEMA P-value',
        gene=schemaDf$Gene.Symbol,score=schemaDf$P.meta))

dim(scoreDf)


# ----------------------------------------------------------------------------------------
# PPI networks to test

networks <- read.table("../data/ASD_RareVariantAnalysis.networks",header=F)$V1
# 29 networks: 13 BAIT, COMBINED, 13 BAIT_noInWeb, COMBINED_noInWeb, InWeb

names <- networks

# read in interactors and non-interactors for each network
intDf <- read.table('../data/ASD_MasterInteractorTable_withInWeb.txt',
	header=T,sep='\t',stringsAsFactors=F)

baits <- sort(unique(intDf$Bait))
baits <- baits[baits != 'COMBINED']

intDf <- subset(intDf, Dataset %in% networks)[,c('Dataset','Gene','IsInteractor')]
dim(intDf)

# network size (int count)
sizeDf <- subset(intDf, IsInteractor)
sizeDf <- data.frame(table(sizeDf$Dataset))
colnames(sizeDf) <- c('Dataset','IntCount')


# ----------------------------------------------------------------------------------------
# repeat anlaysis for pLI, ASD, DD, and SCZ scores

resultsTable <- NULL # enrichment results

for (v in unique(scoreDf$variable)) {

	geneDf <- subset(scoreDf,variable==v)
	
	# KS test direction
	testDir <- NULL
	if (v=='gnomAD pLI') {testDir <- 'less'}
	else {testDir <- 'greater'}
	# 'less': score CDF of x lies BELOW that of y (so x are > y)
	# 'greater': score CDF of x lies ABOVE that of y (so x are < y)
	

	# ------------------------------------------------------------------------------------
	# calculate CONDITIONAL and GLOBAL enrichment for networks
	# (vs. non-interactors and vs. genomic background)
	
	for (i in 1:length(networks)) {

		ints <- subset(intDf,Dataset==networks[i] & IsInteractor)$Gene
		nonInts <- subset(intDf,Dataset==networks[i] & !IsInteractor)$Gene
	
		intScores <- geneDf$score[geneDf$gene %in% ints]
		nonIntScores <- geneDf$score[geneDf$gene %in% nonInts]
		genomeScores <- geneDf$score[!geneDf$gene %in% ints]
		
		# KS test for ints vs. non-ints
		res <- ks.test(intScores,nonIntScores,alternative=testDir)
		
		resultsTable <- rbind(resultsTable, data.frame(Variable=v,Test='Conditional',Set=names[i],
			IntCount=length(intScores),nonIntCount=length(nonIntScores),pDir=testDir,
			Stat=res$statistic,Pvalue=res$p.value))
		
		# KS test for ints vs. genomic background
		res <- ks.test(intScores,genomeScores,alternative=testDir)

		resultsTable <- rbind(resultsTable, data.frame(Variable=v,Test='Global',Set=names[i],
			IntCount=length(intScores),nonIntCount=length(genomeScores),pDir=testDir,
			Stat=res$statistic,Pvalue=res$p.value))

	}	
}


# output enrichment results table
write.table(resultsTable,"../output/ASD_RareVariantEnrichment.txt",
	quote=F,row.names=F,sep="\t")


# ----------------------------------------------------------------------------------------
# plot enrichment results

df <- subset(resultsTable,Set %in% names & Set != 'InWeb')

# set network label and plotting order
df$IntCount <- sizeDf$IntCount[match(df$Set,sizeDf$Dataset)]
df$Bait <- sapply(strsplit(as.character(df$Set),'_'),'[',1) # strip _noInWeb suffix

netOrder <- rev(c('COMBINED','COMBINED_noInWeb',
	sort(names[names!='InWeb' & names!='COMBINED' & names!='COMBINED_noInWeb'])))
netOrder <- unique(paste(df$Bait[match(netOrder,df$Set)],
	' (',sizeDf$IntCount[match(netOrder,sizeDf$Dataset)],')',sep=''))

df$Bait <- paste(df$Bait,' (',df$IntCount,')',sep='')  # add network size
df$Bait <- factor(df$Bait,levels=netOrder)

# set facet order
df$Source <- "This study"
df$Source[grepl("noInWeb",df$Set)] <- "This study - InWeb"
sourceOrder <- c("This study","This study - InWeb")
df$Source <- factor(df$Source, levels=sourceOrder)

df$Test <- factor(df$Test, levels=c('Global','Conditional'))

# set trait order
df$Score[df$Variable=='ASD Q-value'] <- 'ASD'
df$Score[df$Variable=='DD P-value'] <- 'DD'
df$Score[df$Variable=='SCHEMA P-value'] <- 'SCZ'
df$Score[df$Variable=='gnomAD pLI'] <- 'pLI'
df$Score <- factor(df$Score, levels=c('ASD','DD','SCZ','pLI'))


# HEAT MAPS

# nominal and Bonferroni significance (adjusting for 14 networks)
df$Sig <- ifelse(df$Pvalue < 0.05,"*","")
df$Sig[df$Pvalue<0.05/14] <- "**"

df$SigText <- ifelse(df$Pvalue < 0.05,
	paste(formatC(df$Stat,format="f",digits=2),df$Sig,sep=''),"")

pdf("../output/ASD_RareVariantEnrichment_HeatMaps.pdf",height=3,width=8)

# global
ggplot(subset(df,Test=='Global'),aes(x=Score,y=Bait,fill=-log10(Pvalue))) + facet_wrap(~Source,scales='free_y') + 
geom_tile() + geom_text(aes(label=SigText),size=3) +
scale_fill_gradient(name=bquote(-log[10]*"("*italic(.("P"))*")"),low="white",high="red") +
xlab("Score") + ylab("Network") + ggtitle('Global') +
theme_bw() + theme(legend.key.size=unit(0.7,"line"),
	legend.title=element_text(size=9),legend.text=element_text(size=8))

# conditional
ggplot(subset(df,Test=='Conditional'),aes(x=Score,y=Bait,fill=-log10(Pvalue))) + facet_wrap(~Source,scales='free_y') + 
geom_tile() + geom_text(aes(label=SigText),size=3) +
scale_fill_gradient(name=bquote(-log[10]*"("*italic(.("P"))*")"),low="white",high="red") +
xlab("Score") + ylab("Network") + ggtitle('Conditional') + 
theme_bw() + theme(legend.key.size=unit(0.7,"line"),
	legend.title=element_text(size=9),legend.text=element_text(size=8))

dev.off()


# BAR PLOT (conditional, this study)

df$Significant <- "No"
df$Significant[df$Pvalue<0.05] <- "P < 0.05"
df$Significant[df$Pvalue<0.05/14] <- "P < 0.05/14"
sigOrder <- c("P < 0.05/14","P < 0.05","No")
df$Significant <- factor(df$Significant, levels=rev(sigOrder))

df$BarText <- ifelse(df$Pvalue < 0.05,formatC(df$Stat,format="f",digits=2),'')
df$Hjust <- ifelse(grepl('COMBINED',df$Set),1.1,-0.2)
# text alignment: 1=to the left of end of bar, 0=to the right of bar

pdf("../output/ASD_RareVariantEnrichment_BarPlots.pdf",height=3,width=7)

# conditoinal test, this study
ggplot(subset(df,Test=='Conditional' & Source=='This study'),aes(x=Bait,y=-log10(Pvalue),fill=Significant)) +
facet_wrap(~Score,ncol=4,scales='free_x') + geom_col(show.legend=F) +

# significance cutoff lines (behind bars)
geom_hline(yintercept=-log10(0.05),linetype='dashed',size=0.25) +
geom_hline(yintercept=-log10(0.05/14),linetype='dashed',size=0.25) +

# KS statistics (for nominally significant results)
geom_text(aes(label=BarText,hjust=Hjust),size=3) +

scale_fill_manual(values = c('No'='black','P < 0.05'='sienna1','P < 0.05/14'='red')) +
xlab('Network') + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
coord_flip() + theme_bw()

dev.off()

