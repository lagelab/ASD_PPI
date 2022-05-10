##########################################################################################
## Generate GTEx tissue enrichment plots
##
## Author: Kalliopi Tsafou, Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)


# ----------------------------------------------------------------------------------------
# all tissue enrichment
allDf <- read.table('../data/GTEx/ASD_GTExEnrichment_AllRegions.txt',
	header=T,sep="\t",stringsAsFactors=F)
colnames(allDf)[10] <- "Q" # FDR/q-value

# set significance annotation
nTissues <- length(unique(allDf$TissueName)) # 53
allDf$Significant <- ifelse(allDf$p.value < 0.05,'Nominal','None')
allDf$Significant[allDf$p.value<0.05/nTissues] <- 'Bonferroni'

# set tissue plotting order, byt tissue category, then by p-value
allDf <- allDf[order(allDf$TissueCategory,allDf$p.value,decreasing=T),]
allDf$TissueName <- factor(allDf$TissueName, levels=allDf$TissueName)


pdf ('../output/ASD_GTExEnrichment_BarPlots.pdf',height=2.77,width=3)

ggplot(allDf,aes(x=TissueName,y=-log10(p.value),fill=TissueCategory)) + geom_col() +

geom_hline(yintercept=-log10(0.05),linetype='dashed',size=0.25) +
geom_hline(yintercept=-log10(0.05/nTissues),linetype='dashed',size=0.25) +

xlab('Tissue') + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) + 
labs(fill='Tissue category') +
coord_flip() + theme_bw() + theme(axis.title=element_text(size=9),
	axis.text.y=element_blank(),axis.text.x=element_text(size=8),
	legend.text=element_text(size=7),legend.key.size=unit(0.25,'lines'))


# ----------------------------------------------------------------------------------------
# brain tissue enrichment
brainDf <- read.table('../data/GTEx/ASD_GTExEnrichment_BrainRegions.txt',
	header=T,sep="\t",stringsAsFactors=F)
colnames(brainDf)[10] <- "Q" # FDR/q-value

# set significance annotation
nBrains <- length(unique(brainDf$TissueName)) # 13
brainDf$Significant <- ifelse(brainDf$p.value < 0.05,'Nominal','None')
brainDf$Significant[brainDf$p.value<0.05/nBrains] <- 'Bonferroni'

brainDf$SigText <- ifelse(brainDf$p.value < 0.05,brainDf$ip.overlap,'')

# set tissue plotting order, byt tissue category, then by p-value
brainDf$TissueName[brainDf$TissueName=='Cerebellar Hemisphere'] <- 'Cerebellar hemisphere'
brainDf$TissueName[brainDf$TissueName=='Frontal Cortex'] <- 'Frontal cortex'
brainDf <- brainDf[order(brainDf$p.value,decreasing=T),]
brainDf$TissueName <- factor(brainDf$TissueName, levels=brainDf$TissueName)


ggplot(brainDf,aes(x=TissueName,y=-log10(p.value),fill=Significant)) + geom_col() +

geom_hline(yintercept=-log10(0.05),linetype='dashed',size=0.25) +
geom_hline(yintercept=-log10(0.05/nBrains),linetype='dashed',size=0.25) +

geom_text(aes(label=SigText),hjust=-0.1,size=3) +

scale_fill_manual(values = c('None'='black','Nominal'='sienna1','Bonferroni'='red')) +
scale_y_continuous(expand=expansion(mult=c(0,0.2),add=c(0,0))) +

xlab('Brain tissue') + ylab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
coord_flip() + theme_bw() + theme(legend.position='none',
	axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()
