##########################################################################################
## Generate MAGMA common variant enrichment plots for ASD PPI networks
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)


# read in network size
intDf <- read.table('../data/ASD_MasterInteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)

sizeDf <- subset(intDf, IsInteractor)
sizeDf <- data.frame(table(sizeDf$Dataset))
colnames(sizeDf) <- c('Dataset','IntCount')


# read in MAGMA results
df <- read.table('../data/ASD_MagmaResults.txt',header=T,stringsAsFactors=F)
df <- subset(df, (!grepl('_',VARIABLE) & VARIABLE!="InWeb"))


# set network plotting order
df$Network <- df$VARIABLE
baitOrder <- sort(unique(df$Network))
df$Network <- factor(df$Network,
	levels=rev(c("COMBINED",baitOrder[baitOrder!="COMBINED"])))

# add network size to network label
df$IntCount <- sizeDf$IntCount[match(df$Network,sizeDf$Dataset)]

df <- df[order(df$Network),] # factor order already set above
df$Network <- paste(df$Network,' (',df$IntCount,')',sep='')
df$Network <- factor(df$Network,levels=unique(df$Network))


# set trait plotting order
df$TRAIT <- factor(df$TRAIT, levels=c("ASD","ADHD","BIP","MDD","SCZ","height"))


# set significance marker text for nominal for Bonferroni significance
# Bonferroni: adjust for 14 networks
df$Sig <- ifelse(df$P < 0.05,"*","")
df$Sig[df$P<0.05/14] <- "**"

df$SigText <- ifelse(df$P < 0.05,paste(formatC(df$BETA,digits=2),df$Sig,sep=''),'')


# set facet order
df$TEST <- factor(df$TEST,levels=c('Global','Conditional'))


# heat maps with colors base on -log10(P), global vs. conditional as facets
pdf('../output/ASD_MagmaResults_FacetHeatMap.pdf',height=3,width=8)

ggplot(df,aes(x=TRAIT,y=Network,fill=-log10(P))) + facet_wrap(~TEST) + 
geom_tile() + geom_text(aes(label=SigText),size=3) +
scale_fill_gradient(name=bquote(-log[10]*"("*italic(.("P"))*")"),
	low="white", high="red") +
xlab('GWAS trait') + 
theme_bw() + theme(legend.key.size=unit(0.7,"line"),
	axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()

