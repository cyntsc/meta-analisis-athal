##########################################################
###
### Goal: WGCNA Std for A.Thaliana 18 samples
### Method: TOM calculation: adjacency  ( UNSIGNED NTW / Pearson corr)
###
###
### Made by: Cynthia Soto
### Date: December 04, 2020
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by 18 samples of A. Thaliana: ** H E A L T H Y ** & ** I N F E C T E D ** 
### 2) Row RNASeq data comes from SRA-NCBI
###    Data were cleaned, aligned with STAR and quantified with HTSeq-count Tool
### 4) Data were tested with several quality tests before to include in the expression matrix 
###    FastQC / HTSeq-qa / Stats: Central-Measurements & descriptive stat
### 5) Zeros at 100% across all samples were PREVIOUSLY removed  

### Last update: Jan 11, 2021
### Reference: http://pklab.med.harvard.edu/scw2014/WGCNA.html 

##########################################################

#getwd();        ## setwd('../../data/run/cyntsc');
rm(list = ls());

#install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
#source("http://bioconductor.org/biocLite.R") 
#biocLite("impute")

## This library contains all what is nedeed to use herarquical clustering
library(dynamicTreeCut)
library(flashClust)
library(WGCNA);
## The here package allows you to set the top level of your project folder as “here” and to specify where things live relative to that location
library(here);   
here();
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

## Allow multi-treads and setting this in the bash environ
allowWGCNAThreads();
ALLOW_WGCNA_THREADS=12;
## Initial variables 
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

## Load data
here("files", "data", "all_log2_tidy.csv");
athalData <- read.csv(here("data", "all_log2_tidy.csv"), header=TRUE, row.names='Genes', sep='\t')
dim(athalData);
names(athalData);
## Some useless columns and 2 samples were removed from the dataset due their low quality 
athalData = subset(athalData, select = -c(Ss30,Ss30.1)); # 24326 genes (rows) x 19 samples (cols)
dim(athalData);
athalData[1:18];
class(athalData);  #The object is a data.frame

## WGCNA requires genes be given in the columns
gene.names=rownames(athalData);
gene.names;
trans.athalData=t(athalData);

## Focus in a small group of genes to test the function
##n=25000;
##datExpr=trans.athalData[,1:n];

datExpr=trans.athalData
dim(datExpr)
## Pull the names according with the size grp
##SubGeneNames=gene.names[1:n]

## https://rdrr.io/cran/WGCNA/man/pickSoftThreshold.html 
## Choosing a soft-threshold to fit a scale-free topology to the network
powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(
  datExpr,
  dataIsExpr = TRUE,
  powerVector = powers,
  corFnc = cor,                              # cor: Fast calculations of Pearson correlation
  corOptions = list(use = 'p'),
  networkType = "unsigned");                 # "unsigned", "signed", "signed hybrid"

# Plot the results
sizeGrWindow(9, 7)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, unsigned R^2",type="n", main = paste("Scale independence. 18 samples."));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.60,col="red")
#abline(h=0.60,col="blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

abline(h=475,col="red");

softPower = 10;

## AQUI VOY Be sure to change dataframe to matrix type and handle NaN properly (Add.CSC) ############################################

##Error in hclust(d, method, members) : 
##  NA/NaN/Inf in foreign function call (arg 11)

is.matrix(datExpr);
datExpr = data.matrix(datExpr, rownames.force = NA)


####################################################################################################

#calclute the adjacency matrix
adj= adjacency(datExpr,type = "unsigned", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);
## TOM calculation: adjacency..
## ..will use 11 parallel threads.
##  Fraction of slow calculations: 0.000000
## ..connectivity..
## ..matrix multiplication..
## ..normalization..
## ..done.
colnames(TOM) =rownames(TOM) =SubGeneNames
dissTOM=1-TOM

## Module detection
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);
