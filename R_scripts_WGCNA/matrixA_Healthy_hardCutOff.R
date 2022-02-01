##########################################################
###
### Goal: WGCNA Std for A.Thaliana H E A L T H Y
### Made by: Cynthia Soto
### Date: November 09th, 2020
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by 8 samples of A. Thaliana ** HEALTHY** 
### 2) Row RNASeq data comes from SRA-NCBI
###    Data were cleaned, aligned with STAR and quantified with HTSeq-count Tool
### 4) Data were tested with several quality tests before to include in the expression matrix 
###    FastQC / HTSeq-qa / Stats: Central-Measurements & descriptive stat
### 5) Zeros at 100% across all samples were PREVIOUSLY removed  

### Last update: 13 Nov, 2020

##########################################################

getwd();        ## setwd('../../data/run/cyntsc');
#rm(list = ls());

## This library contains all what is nedeed to use herarquical clustering
library(WGCNA);
## The here package allows you to set the top level of your project folder as “here” and to specify where things live relative to that location
library(here);   
here();
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

## Initial variables 
options(stringsAsFactors = FALSE);

## Load data
here("files", "data", "all_log2_dropped20ceros.csv");

athalDataHealthy <- read.csv(here("data", "all_log2_dropped20ceros.csv"), header=TRUE, sep='\t')
#athalDataHealthy;
dim(athalDataHealthy);
names(athalDataHealthy);
athalDataHealthy[1:3]
## Useless columns and control samples are removed from the dataset.
athalDataHealthy = subset(athalDataHealthy, select = -c(X,zero.counter)); 
athalDataHealthy = subset(athalDataHealthy, select =c(1:9));  
dim(athalDataHealthy);   ## 24326 (genes) x 13 (samples)
names(athalDataHealthy);

## remove the auxiliary data 
datExpr0_healthy = as.data.frame(athalDataHealthy[, -c(1:1)]); ## (sample by genes)
## This are the sample names
names(datExpr0_healthy);  ## samples
## This are the gene names
rownames(datExpr0_healthy) = athalDataHealthy$Genes;
rownames(datExpr0_healthy);

#We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0_healthy, verbose = 3);
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts.  If not, we remove the offending genes and samples from the data

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0_healthy)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0_healthy)[!gsg$goodSamples], collapse = ", ")));# Remove the offending genes and samples from the data:
  datExpr0_healthy = datExpr0_healthy[gsg$goodSamples, gsg$goodGenes]	
}


## Build a dendogram to visualize the clusters 
sampleTree = hclust(dist(datExpr0_healthy), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9);

# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#These dataset is 88% from the original size. Zeros across all samples were removed. 
plot(sampleTree, main = "A.thaliana Healthy (8 samples/threshold=17/15.5/9.66)", sub="-", xlab="Genes", ylab = "Height (log2)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

## Set a cutoff line  
abline(h = 17, col = "red");   ## 4 clust (arbitrary threshold selection)
abline(h = 15.5, col = "green");   ## 6 clust (arbitrary threshold selection)
abline(h = 9.66, col="blue");     ## 24 clust / Based on the interquartile grade of the dataset (Q1-Q3)

clust = cutreeStatic(sampleTree, cutHeight = 9.66, minSize = 5)
table(clust);
## with 0 meaning unassigned. The largest cluster is conventionally labeled 1, the next largest 2, etc.



## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# clust 1 contains the genes we want to keep.
keepSamples = (clust==1)
keepSamples; 
datExpr = datExpr0_healthy[keepSamples, ]
nGenes = nrow(datExpr)
nSamples = ncol(datExpr)

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

