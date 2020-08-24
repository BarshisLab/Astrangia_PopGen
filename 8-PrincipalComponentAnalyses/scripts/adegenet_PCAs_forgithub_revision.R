#Script written by Daniel Barshis, adapted by Hannah Aichelman
#Last updated: 06/2020

#This script takes genepop file formats (.gen) and uses principal components analyses to consider population structure

#load required packages
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("viridis")
library("scales")
library("factoextra")

#set working directory
setwd("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/")

#only the host files are shown here for clarity, but just replaced same files for the symbiont analyses to obtain those PCAs
#note that the high outlier and neutral files were also run separately.

datafile<-read.genepop('clone_removed/coral_279_cloneremoved_neutral.filtered1SNPper_genepop.gen', ncode=2)
#datafile<-read.genepop('clone_removed/coral_66_cloneremoved_highoutliers.filtered1SNPper_genepop.gen', ncode=2)


sum(is.na(datafile$tab))
datafile #shows info
YOURdata<-scaleGen(datafile, NA.method='mean')
X<-YOURdata
Y<-as.factor(substring(pop(datafile),1,4)) #change the second number here to change how many letters/numbers you keep in the name
pca1 <- dudi.pca(X,cent=T, scale=T, scannf=F, nf=3)
summary(pca1)

#visualize percent of variance explained by each principal component
fviz_eig(pca1, addlabels=TRUE)

#color symbols, pop names
pdf("HEA_ColorPCA1v2.pdf")
#col <- c("steelblue4","steelblue1", "tomato4", "tomato1") #for coral analysis
col <- c("steelblue4", "tomato4") #for sym analysis

#order is RI_B, RI_W, VA_B, VA_W for coral
#order is RI_B, VA_B for sym

s.class(pca1$li, Y,xax=1,yax=2,
        #sub="RI SNPs, High Outlier Loci",
        #possub = "topleft",
        col=transp(col,.9), 
        axesell=F, 
        cstar=0, 
        cpoint=3,
        grid=FALSE,
        addaxes=TRUE)
add.scatter.eig(pca1$eig[1:3], 3,1,2, posi="topleft")

#title("PCA of HEA_data\naxes 1-2")
dev.off()

