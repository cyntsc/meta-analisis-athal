##########################################################
###
### Goal: WGCNA Std / Signed Ntw with dynamic cut off (pearson) 
###
### Method: TOM calculation: adjacency   ( SIGNED NTW / Pearson corr)
###
###
### Made by: Cynthia Soto 
### Date: March 23, 2021
### Latest update: April 23 2021
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) RNA-Seq libraries are                             ** 17  I N F E C T E D     ** 
###                                                       SAMPLES OF A.THALIANA MATRIX D     
###
### 2) Raw counts were generated with HTSeq-count (python) previously to be aligned with the STAR Tool.
### 3) Raw counts were normalized into TPMs after removing genes with zeros in common thru all the samples (~12% of the gene tarjets)
### 4) After that, resulting array were standarized under 3 criterias: i)remotion of x% of zeros across all samples; ii)TMP<xx; iii)TPM>xxxx
### 5) Resulting array from previous step was transform to Log2(TPM)

# getwd();        ## setwd('../../data/run/cyntsc');
# Clear object lists
rm(list = ls());
# Clear plots
if(!is.null(dev.list())) dev.off()

## This libraries are needed to run this code  
#install.packages(c("WGCNA", dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
#source("http://bioconductor.org/biocLite.R") 
#biocLite("impute")

## import the libraries 
library(tidyverse)
library(dynamicTreeCut)    ## contains all what is nedeed to use herarquical clustering
library(flashClust)
library(WGCNA);
library(lattice, survival, Formula);     ## Let you calculate percentiles and more ...
library(Hmisc);   
library(here);         ## allows you to set the top level of your project folder as “here” and to specify where things live relative to that location   
#here();
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

## Allow multi-treads and setting this in the bash environ
allowWGCNAThreads();
ALLOW_WGCNA_THREADS=30;
## Initial variables 
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

athalData3 <-  read.csv(here("data", "matrix_D_infected.csv"), header=TRUE,  sep='\t')
dim(athalData3);   # 20274  x  18
## Pull the names according with the size grp
SubGeneNames=names(athalData3);

## rearrange columns accordint to traits meta-data (must match)
col_order <- c("Genes", "Bc12",   "Bc12.1", "Bc18",   "Bc18.1", "Bc24",
               "Bc24.1", "Ch22",   "Ch22.1", "Ch22.2", "Ch22.3",
               "Ch40",   "Ch40.1", "Ch40.2", "Ch40.3", "Ss30", "Ss30.1", "Ss30.2")
athalData3 <- athalData3[, col_order]
## You see that genes are listed in a column named "Genes" and samples are in columns
athalData3[1:5,]

## Manipulate athalData3 so it matches the format WGCNA needs
row.names(athalData3) = athalData3$Genes   ## Pull the gene names 
#gene.names=rownames(athalData3)

## Transpose data to set as WGCNA needs
athalData3$Genes = NULL             ## Set genes column as index row
athalData3 = as.data.frame(t(athalData3)) # transpose dataset, now samples are rows and genes are columns
#dim(athalData3) 
athalData3[,1:5]
rownames(athalData3)   ## these are the sample names

## Run this to check if there are gene outliers
gsg = goodSamplesGenes(athalData3, verbose = 3)
gsg$allOK

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK)
   {if (sum(!gsg$goodGenes)>0)
       printFlush(paste("Removing genes:", paste(names(athalData3)[!gsg$goodGenes], collapse= ", ")));
       if (sum(!gsg$goodSamples)>0)
           printFlush(paste("Removing samples:", paste(rownames(athalData3)[!gsg$goodSamples], collapse=", ")))
       athalData3= athalData3[gsg$goodSamples, gsg$goodGenes]
       }

datTraits <- read.csv(here("data", "Traits_Infected_ConditionV2.csv"), header=TRUE,  sep='\t')   
datTraits[1:5,]
datTraits = subset(datTraits, select = -c(accesion, B_cinerea, C_higginsianum, S_sclerotium));

## form a data frame analogous to expression data that will hold the infection traits.
rownames(datTraits) = datTraits$ID
datTraits$ID = NULL                 ## Set ID sample column as index row
table(rownames(datTraits)==rownames(athalData3)) #should return TRUE if datasets align correctly, otherwise your names are out of order

save(athalData3, datTraits, file="Athal_Infected_Log2TPM_SamplesAndTraits.RData")
load("Athal_Infected_Log2TPM_SamplesAndTraits.RData")

###########################################################################################################
#
#                    pickSoftThreshold function offers a analysis of scale free topology for soft-thresholding 
#
##########################################################################################################

## https://rdrr.io/cran/WGCNA/man/pickSoftThreshold.html 
## Choosing a soft-threshold to fit a scale-free topology to the network
#powers = c(c(1:10), seq(from = 12, to=20, by=2));
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
## remember that a signed ntw preserve the natural continuity of the correlation (+ / -), contrary to whats an unsigned ntw does
sft=pickSoftThreshold(athalData3,
  dataIsExpr = TRUE,
  #RsquaredCut = 0.85,                      # desired minimum scale free topology fitting index R^2.
  powerVector = powers,
  corFnc = cor,                             # cor: Fast calculations of Pearson correlation
  corOptions = list(use = 'p'),             # Almost all lists in R internally are Generic Vectors, whereas traditional dotted pair lists (as in LISP) remain available but rarely seen by users (except as formals of functions). 
  verbose =5,
  networkType = "signed");                  # "unsigned", "signed", "signed hybrid"

#warnings()
sft
write.table(sft, file = "results/SFT_corr_Athal_Infected_MatrixD.txt", sep = ",", quote = FALSE, row.names = T)

##Plot the results
sizeGrWindow(5, 5)
par(mfrow = c(1,3));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, Signed R^2 (Log2(TPM))",
     type="n", main = paste("Scale independence for Athal.infected \n (Matrix D)"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red");  

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=1161,col="blue");

# Median connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,6],xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",main = paste("Median connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
abline(h=876 ,col="green");

softPower = 27 # all_infected NO Ss30 samples

# Max correlation factor = 0.8206






