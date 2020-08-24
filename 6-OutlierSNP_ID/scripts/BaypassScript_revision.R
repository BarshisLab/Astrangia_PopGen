# This script was modified by by Hannah Aichelman
# Last updated: 06/2020

# This script was adapted from Sandrine Boissel's wiki page, which was originally taken from the BayPass manual. It runs BayPass (outlier SNP scanning) and analyzes the output. 
# Wiki link here:
# http://192.168.116.226/~barshislab/dokuwiki/doku.php?id=lab_notebooks:sandrine_s_notebook:rerun_snp_filtering_and_baypass

#set working directory and load required packages
setwd("~/Documents/ODU_MS/MolecularWork/RNASeq_Analysis/baypass_2.1/sources")
source("~/Documents/ODU_MS/MolecularWork/RNASeq_Analysis/baypass_2.1/utils/baypass_utils.R")
require(corrplot) ; require(ape)
library(dplyr)

#Note that I have only included the example for one analysis here (coral subset)
#Just changed the file names to analyze the symbiont files and the VA and RI files

#first, upload estimate of omega from baypass
omega=as.matrix(read.table("coral_1808_cloneremoved_revision/coral_1808_cloneremoved_mat_omega.out"))

#Set population names, the order is the same as in the genepop file
pop.names=c("RI_B","VA_W","VA_B","RI_W")

#check dimensions of data frames
dimnames(omega)=list(pop.names,pop.names)

#Compute and visualize the correlation matrix
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
          main=expression("Correlation map based on"~hat(Omega)))

#Visualize the correlation matrix as a hierarchical clustering tree
sites.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(sites.tree,type="p",
      main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

#Estimates of the XtX differentiation measures
sites.snp.res=read.table("coral_1808_cloneremoved_revision/coral_1808_cloneremoved_summary_pi_xtx.out",h=T)

#plot estimates of the XtX differentiation measures
plot(sites.snp.res$M_XtX)

#######################################################
# Generate simulated pseudo-observed data set 
#######################################################
#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution
pi.beta.coef=read.table("coral_1808_cloneremoved_revision/coral_1808_cloneremoved_summary_beta_params.out",h=T)$Mean

#upload the original data to obtain total allele count
Original.data<-geno2YN("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/coral_1808.HWE.final.cloneremoved.recode_baypass.txt")

#Create the POD (pseudo-observed data set)
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=10000,sample.size=Original.data$NN, beta.pi=pi.beta.coef,pi.maf=0,suffix="coral_1808_cloneremoved_POD10kAstrangia") #saves output files with the indicated suffix 

# 1000 SNP simulated out of 10000 
# 2000 SNP simulated out of 10000 
# 3000 SNP simulated out of 10000 
# 4000 SNP simulated out of 10000 
# 5000 SNP simulated out of 10000 
# 6000 SNP simulated out of 10000 
# 7000 SNP simulated out of 10000 
# 8000 SNP simulated out of 10000 
# 9000 SNP simulated out of 10000 
# 10000 SNP simulated out of 10000 

pi.beta.coef
#[1] 5.361749 1.375577  #coral, clone removed

#######################################################
# After POD run through core model
# Sanity Check: Compare POD and original data estimates
# Do this after re-running ./g_baypass script on simulated dataset
#######################################################
#get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table("10KPOD_coral_1808_cloneremoved_revision/10KPOD_coral_1808_cloneremoved_mat_omega.out"))

#plot the estimate and the actual omega
plot(pod.omega,omega) ; abline(a=0,b=1)
fmd.dist(pod.omega,omega) #evaluates distance between two covariance matrices

#[1] 0.4645702  #coral, clone removed

#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table("10KPOD_coral_1808_cloneremoved_revision/10KPOD_coral_1808_cloneremoved_summary_beta_params.out",h=T)$Mean

plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)

#######################################################
# XtX calibration
#######################################################
#get the pod XtX

pod.xtx=read.table("10KPOD_coral_1808_cloneremoved_revision/10KPOD_coral_1808_cloneremoved_summary_pi_xtx.out", h=T)$M_XtX

#compute the 1% threshold
pod.thresh=quantile(pod.xtx,probs=0.99)

#add the threshold to the actual XtX plot
plot(sites.snp.res$M_XtX)
abline(h=pod.thresh,lty=2)


#############here is where the BayPass manual example stops under the core model mode##############
#adapted code below from Sandrine Boissel's wiki
#prevsites.snp.res=read.table("Carolyn_results/pmecore_summary_pi_xtx.out",h=T)

sites.snp.res=read.table("coral_1808_cloneremoved_revision/coral_1808_cloneremoved_summary_pi_xtx.out",h=T)

plot(sites.snp.res$M_XtX)

pod.xtx=read.table("10KPOD_coral_1808_cloneremoved_revision/10KPOD_coral_1808_cloneremoved_summary_pi_xtx.out", h=T)$M_XtX

#compute the 0.5% threshold
pod.highthresh=quantile(pod.xtx,probs=0.95) #"high" outlier SNPs
pod.lowthresh=quantile(pod.xtx,probs=0.05) #low outlier SNPs
pod.higherthresh=quantile(pod.xtx,probs=0.99) #"higher" outlier SNPs
pod.lowerthresh=quantile(pod.xtx,probs=0.01) #"lower" outlier SNPs

#add the threshold to the actual XtX plot
sites.snp.res$color="black"
sites.snp.res$color[sites.snp.res$M_XtX>pod.highthresh]="red"
sites.snp.res$color[sites.snp.res$M_XtX<pod.lowthresh]="blue"
sites.snp.res$SNP="neutral"
sites.snp.res$SNP[sites.snp.res$M_XtX>pod.highthresh]="high"
sites.snp.res$SNP[sites.snp.res$M_XtX<pod.lowthresh]="low"
sites.snp.res$SNP[sites.snp.res$M_XtX>pod.higherthresh]="higher"
sites.snp.res$SNP[sites.snp.res$M_XtX<pod.lowerthresh]="lower"
plot(sites.snp.res$M_XtX, col=sites.snp.res$color)
abline(h=c(pod.highthresh,pod.higherthresh,pod.lowthresh,pod.lowerthresh),lty=c(1,2,1,2),col=c("red","red","blue","blue"))

#write out a csv of the outlier SNPs
write.csv(sites.snp.res, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/RI_revision_outlier_SNPs.csv")

#The rest of the script creates files that include either all SNPs, only high outlier SNPs, or only low outlier SNPs. 
#These files were used to subset the vcf files

#Outlier files:
#The contig names here come from the genepop file, the contigs are in the same order as the MRK column output above 
coralSNP<-read.csv("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/coral_1808_cloneremoved_SNPslabeled.csv")

#coralSNP<-read.csv("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/VA_revision_outlier_SNPsLabeled.csv")

coralPopOut <- coralSNP %>%
  select(Contig, MRK, SNP) %>%
  filter(SNP=="high"|SNP=="higher"|SNP=="low"|SNP=="lower")%>%
  select(Contig)
names(coralPopOut) <- NULL
head(coralPopOut)
write.table(coralPopOut, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/VA_AllOutliers.txt", row.names=FALSE, quote=FALSE)

coralHighOut <- coralSNP %>%
  dplyr::select(Contig, MRK, SNP) %>%
  filter(SNP=="high"|SNP=="higher")%>%
  dplyr::select(Contig)
names(coralHighOut) <- NULL
head(coralHighOut)
write.table(coralHighOut, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/VA_HighOutliers.txt", row.names=FALSE, quote=FALSE)

coralLowOut <- coralSNP %>%
  dplyr::select(Contig, MRK, SNP) %>%
  filter(SNP=="low"|SNP=="lower")%>%
  dplyr::select(Contig)
names(coralLowOut) <- NULL
head(coralLowOut)
write.table(coralLowOut, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/VA_LowOutliers.txt", row.names=FALSE, quote=FALSE)

coralNeutral <- coralSNP %>%
  dplyr::select(Contig, MRK, SNP) %>%
  filter(SNP=="neutral")%>%
  dplyr::select(Contig)
names(coralNeutral) <- NULL
head(coralNeutral)
write.table(coralNeutral, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/VA_Neutral.txt", row.names=FALSE, quote=FALSE)


#Look for overlap in VA and RI high outliers
RI_highout <- read.delim("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/RI_HighOutliers.txt")
VA_highout <- read.delim("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/VA_HighOutliers.txt")

dim(RI_highout) #47 high outliers
dim(VA_highout) #49 high outliers

match<-match_df(RI_highout, VA_highout, on="Contig")
match2<-match_df(VA_highout, RI_highout, on="Contig")
#returns the matching rows from the dataset listed first. We get the same contigs in both lists, just different numbers of both
write.table(match, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/VA_RI/VA_RI_matchinghighoutliers.txt", row.names=FALSE, quote=FALSE)


#Outlier files for brown individual symbiont subset analyses:
SymSNP<-read.csv("/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/sym_59_SNPslabeled.csv")

SymOut <- SymSNP %>%
  dplyr::select(Contig, MRK, SNP) %>%
  filter(SNP=="high"|SNP=="higher"|SNP=="low"|SNP=="lower")%>%
  dplyr::select(Contig)
names(SymOut) <- NULL
head(SymOut)
write.table(SymOut, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/sym_AllOutliers.txt", row.names=FALSE, quote=FALSE)

SymHighOut <- SymSNP %>%
  dplyr::select(Contig, MRK, SNP) %>%
  filter(SNP=="high"|SNP=="higher")%>%
  dplyr::select(Contig)
names(SymHighOut) <- NULL
head(SymHighOut)
write.table(SymHighOut, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/sym_HighOutliers.txt", row.names=FALSE, quote=FALSE)

SymLowOut <- SymSNP %>%
  dplyr::select(Contig, MRK, SNP) %>%
  filter(SNP=="low"|SNP=="lower")%>%
  dplyr::select(Contig)
names(SymLowOut) <- NULL
head(SymLowOut)
write.table(SymLowOut, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/sym_LowOutliers.txt", row.names=FALSE, quote=FALSE)

SymNeutral <- SymSNP %>%
  dplyr::select(Contig, MRK, SNP) %>%
  filter(SNP=="neutral")%>%
  dplyr::select(Contig)
names(SymNeutral) <- NULL
head(SymNeutral)
write.table(SymNeutral, file="/Users/hannahaichelman/Documents/Manuscripts/AstrangiaSNPs/ForGithub/5-OutlierSNP_ID/revision_files/sym_Neutral.txt", row.names=FALSE, quote=FALSE)
