library(vegan)
library(ggplot2)
library(grid)
library(plyr)

GeneraK <- read.csv("GeneraKraken.csv",header=TRUE,row.names=1)
colnames(GeneraK) <- gsub("_Sub","",colnames(GeneraK))
GeneraK <- t(GeneraK)
GeneraP <- GeneraK/rowSums(GeneraK)
Meta <- read.csv("../data/metaFP1B.csv",header=TRUE,row.names=1)

GeneraP <- GeneraP[rownames(Meta),]
GeneraP.nmds <- metaMDS(GeneraP)

adonis(GeneraP ~ Meta$Reactor)

adonis(GeneraP ~ Meta$Time)


nmds_df<-scores(GeneraP.nmds,display=c("sites"))

nmds_df<-data.frame(nmds_df)

meta_nmds.df <- data.frame(x=nmds_df$NMDS1,y=nmds_df$NMDS2,Day=Meta$Day,Reactor=as.factor(Meta$Reactor))

reactor_colours = c(R1 = "red", R2 = "green", R3 ="blue")


p<-ggplot(data=meta_nmds.df,aes(x,y,colour=Reactor)) + geom_point(alpha=0.5) + scale_colour_manual(values=reactor_colours)

pdf("GeneraKNMDS.pdf")
plot(p + geom_path(arrow=arrow(length=unit(0.3,"cm")),alpha=0.5) + theme_bw())
dev.off()
