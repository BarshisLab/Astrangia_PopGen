#This script reads in vcf files and uses the package poppr to detect clones

#Written by Hannah Aichelman
#Last updated: July 2020

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ade4)
library(adegenet)
library(mmod)
library(strataG)
library(hierfstat)
library(ape)
library(StAMPP)
library(pegas)
library(poppr)
library(Rcpp)
library(factoextra)
library(reshape2)
library(vegan)
library(stringr)
library(dartR)
library(gplots)
library(RColorBrewer)

setwd("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/Revision_PeerJ/Clone_Analysis")

#### using new vcf files after updating SNP filtering during manuscript revision (7/13/2020)
load("Hannah_new_snps_genind_genlight.RData")
#this .Rdata contains two objects, a genind and a genlight. 
#the lines of code to produce these files in the .Rdata are:
#hannah_coral_vcf<-read.vcfR("coral_1808.HWE.final.recode.vcf")
#hannah_coral_genind<-vcfR2genind(hannah_coral_vcf)
#hannah_coral_genlight<-vcfR2genlight(hannah_coral_genlight)


#Dissimilarity Distance
#Create a diss.dist dissimilarity distance matrix to determine a threshhold based on clones, https://www.rdocumentation.org/packages/poppr/versions/2.8.5/topics/diss.dist
#ignores missing data and counts the shared genotypes 
distgenDISS <- diss.dist(hannah_coral_genind, percent = FALSE, mat = FALSE) 
#make percent = TRUE to get Prevosti distance
#By default, diss.dist() returns a distance reflecting the number of allelic differences between two individuals (Hamming's distance).

distgenDISS2 <- as.matrix(distgenDISS)
write.csv(distgenDISS,"distgenDISS.csv")

#make a heat map of genetic/allelic distances
heatmap.2(distgenDISS2, Rowv = TRUE, Colv = TRUE, dendrogram = c("both"), col=brewer.pal(9, "Blues"),trace = c("none"),cexRow = 0.4, cexCol = 0.4, density.info = c("none"), key = TRUE)

#Cluster with hclust and plot
clust_tree<-hclust(distgenDISS, "ave")
plot(clust_tree, cex=0.4)+abline(h=50)


#make a heat map of Manhattan distance - just checking to see if other distance methods give us the same answer
Man <- vegdist(hannah_coral_genind, method = "manhattan", na.rm = TRUE)
Man.matrix <- as.matrix(Man)
#View(Man.matrix)
heatmap.2(Man.matrix, Rowv = TRUE, Colv = TRUE, dendrogram = c("both"), col=brewer.pal(9, "Blues"),trace = c("none"),cexRow = 0.4, cexCol = 0.4, density.info = c("none"), key = TRUE)

#Dissimilarity Distance - Prevosti
distgenDISS_P <- diss.dist(hannah_coral_genind, percent = TRUE, mat = FALSE) 

#Cluster with hclust and plot
clust_tree_p<-hclust(distgenDISS_P, "ave")
plot(clust_tree_p, cex=0.4)+abline(h=50)

## all analyses show RI_B_06_merged and RI_B_04_merged as much more similar than the rest of the samples
## therefore calling these putative clones, and removed (randomly) RI_B_06_merged

#### Stopped here for what is presented in the manuscript, but other ways of looking at clones presented below

#plot genlight object
glPlot(hannah_coral_genlight)

#convert to genind in order to plot cluster dendrogram again
genind<-gl2gi(hannah_coral_genlight)

distgenDISS <- diss.dist(genind, percent = FALSE, mat = FALSE)

clust_tree<-hclust(distgenDISS, "ave")
plot(clust_tree, cex=0.4)+abline(h=23)

#### with new symbiont vcf files after updating SNP filtering (7/13/2020)
load("Hannah_new_sym_snps_genind_genlight.RData")
#this .Rdata contains two objects, a genind and a genlight. 
#the lines of code to produce these files in the .Rdata are:
#hannah_coral_vcf<-read.vcfR("coral_1808.HWE.final.recode.vcf")
#hannah_coral_genind<-vcfR2genind(hannah_coral_vcf)
#hannah_coral_genlight<-vcfR2genlight(hannah_coral_genlight)


#start with the genind object
distgenDISS_sym <- diss.dist(hannah_sym_genind, percent = FALSE, mat = FALSE)

clust_tree<-hclust(distgenDISS_sym, "ave")
plot(clust_tree, cex=0.4)+abline(h=23)
#just shows the VA/RI difference, but too few SNPs to actually call these "clones"

#plot genlight object
glPlot(hannah_sym_genlight)

#convert to genind in order to plot cluster dendrogram again
genind<-gl2gi(hannah_sym_genlight)

distgenDISS <- diss.dist(genind, percent = FALSE, mat = FALSE)

clust_tree<-hclust(distgenDISS, "ave")
plot(clust_tree, cex=0.4)+abline(h=23)
