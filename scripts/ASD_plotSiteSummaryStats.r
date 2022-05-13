##########################################################################################
## Compare and plot summary statistics of IP-MS datasets subsetted by MS site
##
## Author: Yu-Han Hsu, Joshua Ching
##########################################################################################

rm(list=ls())

library(reshape2)
library(ggpubr)


# ----------------------------------------------------------------------------------------
# compare summary stats for the 26 IP-MS datasets that passed QC, subsetted by MS site

df <- read.table('../data/ASD_IpSummaryStats.txt',header=T,sep='\t',stringsAsFactors=F)
df$InWebP <- -log10(df$InWebP)
df$Bait_FDR <- -log10(df$Bait_FDR)
df$Med_FDR <- -log10(df$Med_FDR)
df$Max_FDR <- -log10(df$Max_FDR)

df <- melt(df,id.vars=c('IP','Site','Labeling'))
df <- subset(df, variable %in% c('RepCorr','NumProteins','NumInts','InWebP',
	'Bait_logFC','Bait_FDR','Med_logFC','Med_FDR','Max_logFC','Max_FDR'))

df$variable <- as.character(df$variable)
df$variable[df$variable=='RepCorr'] <- 'Replicate correlation'
df$variable[df$variable=='NumProteins'] <- '# of proteins'
df$variable[df$variable=='NumInts'] <- '# of interactors'
df$variable[df$variable=='InWebP'] <- '-log10 InWeb P'
df$variable[df$variable=='Bait_logFC'] <- 'Bait log2 FC'
df$variable[df$variable=='Bait_FDR'] <- 'Bait -log10 FDR'
df$variable[df$variable=='Med_logFC'] <- 'Median log2 FC'
df$variable[df$variable=='Med_FDR'] <- 'Median -log10 FDR'
df$variable[df$variable=='Max_logFC'] <- 'Max log2 FC'
df$variable[df$variable=='Max_FDR'] <- 'Max -log10 FDR'

df$Site <- factor(df$Site, levels=c('Broad','WH','CNCR'))
df$variable <- factor(df$variable, levels=as.character(unique(df$variable)))


# set Wilcoxon test parameters
varPairs <- combn(levels(df$Site),2)
testList <- NULL
for (col in 1:ncol(varPairs)) { testList[[col]] <- varPairs[,col] }
pCutoffs <- c(0,0.05/length(testList),0.05,2) # set upperbound to 2 (to avoid "NS" label when p=1)
pSymbols <- c('**','*','ns')


# box plots with Wilcoxon p-value annotations
pdf("../output/ASD_SiteSummaryStats_BoxPlots.pdf",height=3,width=8)

ggboxplot(df,x='Site',y='value',color='Site',add='jitter',add.params=list(size=0.7)) +
facet_wrap(~variable,ncol=5,scales='free_y') +

# pairwise Wilcoxon tests
stat_compare_means(comparisons=testList,method='wilcox.test',
symnum.args=list(cutpoints=pCutoffs,symbols=pSymbols),hide.ns=T,size=3,vjust=0.6) +

scale_color_brewer(palette='Set1') +
xlab('MS site') + ylab('Statistic') +
theme_bw() + theme(legend.position='none',axis.text=element_text(size=8))

dev.off()


# ----------------------------------------------------------------------------------------
# compare stats for IP pairs of the same bait, subsetted by same vs different MS site

df <- read.table('../data/ASD_IpComparisonStats.txt',header=T,sep='\t',stringsAsFactors=F)

# MS site of IP1 and IP2
df$Site1 <- unlist(lapply(strsplit(df$IP1,'_'),'[[',2))
df$Site1 <-  unlist(lapply(strsplit(df$Site1,'[0-9]'),'[[',1))
df$Site2 <- unlist(lapply(strsplit(df$IP2,'_'),'[[',2))
df$Site2 <-  unlist(lapply(strsplit(df$Site2,'[0-9]'),'[[',1))

# annotate if IP pair from same or different sites
df$Sites[df$Site1==df$Site2] <- 'Same'
df$Sites[df$Site1!=df$Site2] <- 'Different'
df$Sites <- factor(df$Sites,levels=c('Same','Different'))


# set Wilcoxon test parameters
varPairs <- combn(levels(df$Sites),2)
testList <- NULL
for (col in 1:ncol(varPairs)) { testList[[col]] <- varPairs[,col] }
pCutoffs <- c(0,0.05,2) # set upperbound to 2 (to avoid "NS" label when p=1)
pSymbols <- c('*','ns')

# data frame containing stats to plot
plotDf <- rbind(
        data.frame(Stat='Protein overlap %',
                Sites=c(df$Sites,df$Sites),Value=c(df$Population.IP1,df$Population.IP2)*100),

        data.frame(Stat='log2 FC correlation',Sites=df$Sites,Value=df$logFC_Correlation),

        data.frame(Stat='Interactor overlap %',
                Sites=c(df$Sites,df$Sites),Value=c(df$Overlap.IP1,df$Overlap.IP2)*100),

        data.frame(Stat='-log10 interactor overlap P',Sites=df$Sites,Value=-log10(df$P)))

plotDf$Stat <- factor(plotDf$Stat,levels=unique(plotDf$Stat))


# box plots with Wilcoxon p-value annotations
pdf("../output/ASD_IpPairsBySite_BoxPlots.pdf",height=2,width=8)

ggboxplot(plotDf,x='Sites',y='Value',color='Sites',add='jitter',add.params=list(size=0.7)) +
facet_wrap(~Stat,ncol=4,scales='free_y') +

# pairwise Wilcoxon tests
stat_compare_means(comparisons=testList,method='wilcox.test',
symnum.args=list(cutpoints=pCutoffs,symbols=pSymbols),hide.ns=T,size=3,vjust=0.6) +

scale_color_brewer(palette='Set1') +
xlab('MS sites') + ylab('Statistic') +
theme_bw() + theme(legend.position='none',axis.text=element_text(size=8))

dev.off()
