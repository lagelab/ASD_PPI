##########################################################################################
## Generate network plot and pie charts for ASD PPI data
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(igraph)
library(qgraph)
library(ggplot2)
library(RColorBrewer)
library(genoppi)


# ----------------------------------------------------------------------------------------
# interactor table

intTable <- read.table("../data/ASD_MasterInteractorTable.txt",
	header=T,sep="\t",stringsAsFactors=F)
intTable <- subset(intTable,!grepl("_",Dataset) & Dataset!="COMBINED" &
	IsInteractor)[,c("Bait","Gene")]
colnames(intTable) <- c("Bait","GeneName")


# InWeb interactors of each bait
uniqBaits <- unique(intTable$Bait)
inwebDf <- NULL
for (bait in uniqBaits) {
	tempDf <- genoppi::get_inweb_list(bait)
	if (!is.null(tempDf)) { # if bait in InWeb
		inwebInts <- subset(tempDf,significant==T)$gene
		inwebDf <- rbind(inwebDf,data.frame(Bait=bait,gene=inwebInts))
	}
}

# inInWeb boolean vector (for annotating network edge)
intTable$inInWeb <- NA
for (i in 1:nrow(intTable)) {
	intTable$inInWeb[i] <- intTable$GeneName[i] %in%
		inwebDf$gene[inwebDf$Bait==intTable$Bait[i]]
}


# order columns and rows in intTable to plot bait nodes and InWeb edges last/on top
intTable <- intTable[order(intTable$inInWeb),]

# create network
intNetwork <- graph_from_data_frame(intTable[,c('GeneName','Bait')], directed=F)


# set vertex attributes (for node color, size, label)
V(intNetwork)$isBait <- V(intNetwork)$name %in% unique(intTable$Bait)

maxNumBaits <- max(table(intTable$GeneName))
nodeColors <- brewer.pal(maxNumBaits+1,"Purples")[2:(maxNumBaits+1)]

V(intNetwork)$size <- NA
V(intNetwork)$color <- NA
for (x in V(intNetwork)$name) {
	if (x %in% unique(intTable$Bait)) {
		V(intNetwork)$size[V(intNetwork)$name==x] <- maxNumBaits
		V(intNetwork)$color[V(intNetwork)$name==x] <- "red"
	} else {
		V(intNetwork)$size[V(intNetwork)$name==x] <- sum(intTable$GeneName==x)
		V(intNetwork)$color[V(intNetwork)$name==x] <-
			nodeColors[sum(intTable$GeneName==x)]
	}
}

V(intNetwork)$label <- NA
V(intNetwork)$label[V(intNetwork)$size>=2] <- V(intNetwork)$name[V(intNetwork)$size>=2]


# set edge attributes (for edge color)
E(intNetwork)$inInWeb <- intTable$inInWeb


# network plot layout
e <- get.edgelist(intNetwork,names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(intNetwork),niter=5000,
	area=(vcount(intNetwork)^1.7),repulse.rad=(vcount(intNetwork)^2.2))
			
			
# generate network plot
pdf("../output/ASD_NetworkPlot.pdf",width=7,height=7)

# with gene name labels
plot(intNetwork, edge.width=0.5,
	vertex.color=V(intNetwork)$color,
	vertex.size=V(intNetwork)$size+1, vertex.frame.color="grey30",
	vertex.label=V(intNetwork)$label, vertex.label.color="black",
	vertex.label.cex=0.08,
	edge.color=c("grey","blue")[1+(E(intNetwork)$inInWeb)],
	layout=l)

# without gene name labels
plot(intNetwork, edge.width=0.5,
	vertex.color=V(intNetwork)$color,
	vertex.size=V(intNetwork)$size+1, vertex.frame.color="grey30",
	vertex.label=NA, vertex.label.cex=0.08,
	edge.color=c("grey","blue")[1+(E(intNetwork)$inInWeb)],
	layout=l)

dev.off()


# ----------------------------------------------------------------------------------------
# generate pie charts

# (1) number of index proteins linked to each interactor
genePieTable <- data.frame(table(table(intTable$GeneName)))
colnames(genePieTable) <- c("NumBaits","GeneCount")
genePieTable$Percent <- genePieTable$GeneCount/sum(genePieTable$GeneCount)*100

print(genePieTable)


# (2) interactions in InWeb vs. not in InWeb
inwebTable <- data.frame(table(intTable$inInWeb))
inwebTable$Percent <- inwebTable$Freq/sum(inwebTable$Freq)*100

print(inwebTable)


pdf("../output/ASD_NetowrkPieCharts.pdf",width=3,height=3)

ggplot(genePieTable, aes(x="", y=GeneCount, fill=NumBaits)) +
geom_bar(width=1, stat="identity", color="black") +
coord_polar("y") +
scale_fill_manual(name="Index\nproteins", values=nodeColors) +
theme_void()

ggplot(inwebTable, aes(x="", y=Freq, fill=Var1)) +
geom_bar(width=1, stat="identity", color="black") +
coord_polar("y") +
scale_fill_manual(name="InWeb",labels=c('No','Yes'),values=c("grey","blue")) +
theme_void()

dev.off()

