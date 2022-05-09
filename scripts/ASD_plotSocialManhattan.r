##########################################################################################
## Generate social Manhattan plot using Satterstrom 2020 data
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggrepel)


# ----------------------------------------------------------------------------------------
# read in Satterstrom data
geneDf <- read.table("../data/Satterstrom2020_TableS2_Autosomal.txt",header=T,sep="\t")
geneDf <- geneDf[,c("hugoGene","chr","start_hg38","qval_dnccPTV")]
colnames(geneDf) <- c("Gene","Chr","Pos","Qvalue")
geneDf$Chr <- as.numeric(gsub("chr","",geneDf$Chr))
geneDf[geneDf$Qvalue==0,"Qvalue"] <- 1e-16 # set minimum q-value to avoid Inf in plot 
geneDf <- geneDf[order(geneDf$Chr,geneDf$Pos),]
dim(geneDf)

nChr <- length(unique(geneDf$Chr))
geneDf$Coord <- NA
s <- 0
nbp <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
	138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
	83257441,80373285,58617616,64444167,46709983,50818468)
	# GRCh38.p13 chr 1-22 lengths from: https://www.ncbi.nlm.nih.gov/grc/human/data

chrTicks <- NULL
for (i in 1:22) {
	geneDf[geneDf$Chr==i,"Coord"] <- geneDf[geneDf$Chr==i,"Pos"] + s
	chrTicks <- c(chrTicks,floor(nbp[i]/2) + s)
	s <- s + nbp[i]
}


# ----------------------------------------------------------------------------------------
# read in interactor list
intDf <- read.table("../data/ASD_MasterInteractorTable.txt",header=T,sep="\t")
intDf <- subset(intDf, !grepl("_",Dataset) & Dataset!="COMBINED" & 
	IsInteractor)[,c("Bait","Gene")]
dim(intDf)


# data frame for genes that have qval <= 0.25 and are baits or interactors in ASD IPs
plotDf <- merge(geneDf,intDf,by.x="Gene",by.y="Gene")
plotDf$Bait_Qvalue <- NA
plotDf$Bait_Coord <- NA
baitList <- unique(plotDf$Bait)
for (bait in baitList) {
	# add bait coordinates (for plotting social links between genes and baits)
	plotDf[plotDf$Bait==bait,"Bait_Qvalue"] <- geneDf$Qvalue[geneDf$Gene==bait]
	plotDf[plotDf$Bait==bait,"Bait_Coord"] <- geneDf$Coord[geneDf$Gene==bait]
}
for (bait in baitList) {
	baitDf <- data.frame(geneDf[geneDf$Gene==bait,c("Gene","Chr","Pos","Qvalue","Coord")],
		Bait=NA,Bait_Qvalue=NA,Bait_Coord=NA)
	plotDf <- rbind(plotDf,baitDf,row.names=NULL)
}
dim(plotDf)

# significance group (for node color)
plotDf <- plotDf[plotDf$Qvalue<=0.25,]
plotDf$Significance <- NA
plotDf[plotDf$Qvalue <= 0.25,"Significance"] <- "FDR<=0.25"
plotDf[plotDf$Qvalue <= 0.1,"Significance"] <- "FDR<=0.1"
plotDf[plotDf$Gene %in% baitList,"Significance"] <- "BAIT"

# other plotting parameters
ylim <- abs(floor(log10(min(plotDf$Qvalue)))) + 0.5 

# unique genes that show up in plotDf (which contains unique PAIRS of genes)
uniqDf <- unique(plotDf[,c("Gene","Coord","Qvalue","Significance")])


# ----------------------------------------------------------------------------------------
# generate social Manhttan plot

# edge color based on validation data
wbDf <- read.table('../data/ASD_ValidationSummary.txt',
	header=T,sep='\t',stringsAsFactors=F,row.names=1)
wbDf <- subset(wbDf,Validated)

plotDf$Validated <- FALSE
for (i in 1:nrow(wbDf)) {
	plotDf$Validated[which(plotDf$Gene==wbDf$Interactor[i] &
		plotDf$Bait==wbDf$Bait[i])] <- TRUE
}
plotDf <- plotDf[order(plotDf$Validated),]


# node size based on # of interactions
sizeDf <- data.frame(table(plotDf$Gene))
maxNumBaits <- max(sizeDf$Freq)
uniqDf$NodeSize <- sizeDf$Freq[match(uniqDf$Gene,sizeDf$Var1)]
uniqDf$NodeSize[uniqDf$Gene %in% baitList] <- maxNumBaits


set.seed(1234) # fix ggrepel label coordinates

pdf("../output/ASD_SocialManhattanPlot.pdf",width=7.5,height=3.5)

ggplot(uniqDf,mapping=aes(x=Coord,y=-log10(Qvalue),size=NodeSize,
	color=as.factor(Significance))) + 

# EWS (Q-value for DSCAM) and FDR <= 0.1 line
geom_hline(yintercept=-log10(0.00013466435347903),
	color="black",linetype="dashed",size=0.25) +
geom_hline(yintercept=-log10(0.1),color="black",linetype="dashed",size=0.25) +

# social links
geom_segment(plotDf,mapping=aes(x=Coord,y=-log10(Qvalue),
	xend=Bait_Coord,yend=-log10(Bait_Qvalue)),size=0.3,alpha=0.7,
	color=ifelse(plotDf$Validated,'darkorange','grey')) +

# data point for each gene
geom_point(show.legend=F) +

# gene name label
geom_text_repel(mapping=aes(label=Gene,color=as.factor(Significance),size=NodeSize),
	fontface='bold',segment.size=0.1,segment.alpha=0.5,max.overlaps=100,show.legend=F) +

scale_size(range=c(1.5,2.5)) +
scale_color_manual(name="ASD gene groups",labels=c("Index gene","FDR<=0.1","FDR<=0.25"),
	values=c("red","blue","darkcyan")) +
scale_x_continuous(limits=c(0,max(geneDf$Coord)),label=as.character(1:22),breaks=chrTicks) +
scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) +

#ggtitle("Social Manhattan Plot (Satterstrom et al. 2020)") +
xlab("Chromosomal Position") + ylab(bquote(-log[10]*"("*italic(.("Q"))*"-value)")) +
theme_classic() +
theme(axis.title=element_text(size=9),axis.text=element_text(size=8))

dev.off()

