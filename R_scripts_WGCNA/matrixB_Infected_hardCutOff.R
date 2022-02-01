##########################################################
###
### Goal: WGCNA Std for A.Thaliana I N F E C T E D
### Made by: Cynthia Soto
### Date: November 09th, 2020
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by 12 samples of A. Thaliana I N F E C T E D
###                              ###  ###   ### ##### ##### ##### ###### ##### ## ##
###                              ###  ####  ### ##    ##    ##      ##   ##    ##  ##
###                              ###  ### ###   ####  ####  ##      ##   ###   ##   ##
###                              ###  ###  #### ##    ##    ##      ##   ##    ##   ##
###                              ###  ###   ### ##    ##### ####    ##   ##### ## ##
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

athalDataInfec <- read.csv(here("data", "all_log2_dropped20ceros.csv"), header=TRUE, sep='\t')
#athalDataInfec;
dim(athalDataInfec);
names(athalDataInfec);
athalDataInfec[1:3]
## Useless columns and control samples are removed from the dataset.
athalDataInfec = subset(athalDataInfec, select = -c(X,zero.counter)); 
## Control samples are: SRR3383640,SRR3383641, SRR3383782, SRR3383783, SRR3383821,SRR3383822,SRR6283144,SRR6283145
athalDataInfec = subset(athalDataInfec, select = -c(SRR3383640,SRR3383641, SRR3383782, SRR3383783, SRR3383821,SRR3383822,SRR6283144,SRR6283145)); 
## Some useless columns and 2 samples were removed from the dataset due their low quality
#athalDataInfec = subset(athalDataInfec, select = -c(X,zero.counter,SRR6283147,SRR6283148,zero.counter)); # 24326 genes (rows) x 19 samples (cols)
dim(athalDataInfec);   ## 24326 (genes) x 13 (samples)
names(athalDataInfec);

## remove the auxiliary data and transpose the expression data for further analysis.
datExpr0_Infect = as.data.frame(athalDataInfec[, -c(1:1)]); ## (sample by genes)
#datExpr0_Infect[1:10];

## This are the sample names
names(datExpr0_Infect);  ## samples
## This are the gene names
rownames(datExpr0_Infect) = athalDataInfec$Genes;
rownames(datExpr0_Infect);

#We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0_Infect, verbose = 3);
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts.  If not, we remove the offending genes and samples from the data

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0_Infect)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0_Infect)[!gsg$goodSamples], collapse = ", ")));# Remove the offending genes and samples from the data:
  datExpr0_Infect = datExpr0_Infect[gsg$goodSamples, gsg$goodGenes]	
}


## Build a dendogram to visualize the clusters 
sampleTree = hclust(dist(datExpr0_Infect), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9);

# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#These dataset is 88% from the original size. Zeros across all samples were removed. 
plot(sampleTree, main = "A.thaliana Infected (12 samples/threshold=17/15.5/9.39)", sub="", xlab="Genes", ylab = "Height (log2)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

## Set a cutoff line  
abline(h = 17, col = "red");   ## 4 clust (arbitrary threshold selection)
abline(h = 15.5, col = "green");   ## 6 clust (arbitrary threshold selection)
abline(h = 9.39, col="blue");     ## 24 clust / Based on the interquartile grade of the dataset (Q1-Q3)

clust = cutreeStatic(sampleTree, cutHeight = 9.39, minSize = 5)
# clust = cutreeStatic(sampleTree, cutHeight = 9.24, minSize = 5)
## with 0 meaning unassigned. The largest cluster is conventionally labeled 1, the next largest 2, etc.
table(clust);

# clust 1 contains the genes we want to keep.
keepSamples = (clust==1)
keepSamples; 
datExpr = datExpr0_Infect[keepSamples, ]
nGenes = nrow(datExpr)
nSamples = ncol(datExpr)


### /////////////////////////////////////////////////////////////////////////////////////////////////////////
###                          Data assay 2 on   
###                          Dataset1: Athal infected 
###                          Sclerotinium samples will be removed due low % alignment and atípical KDE curve (not bimodal)
### /////////////////////////////////////////////////////////////////////////////////////////////////////////


athalDataInfec = subset(athalDataInfec, select = -c(SRR6283147,SRR6283148)); 
dim(athalDataInfec);   ## 24326 (genes) x 13 (samples)
#sample names
names(athalDataInfec);

## remove the auxiliary data and transpose the expression data for further analysis.
datExpr0_Infect = as.data.frame(athalDataInfec[, -c(1:1)]); ## (sample by genes)
datExpr0_Infect[1:5];

# This are the sample names
names(datExpr0_Infect);  ## samples
# This are the gene names
rownames(datExpr0_Infect) = athalDataInfec$Genes;
rownames(datExpr0_Infect);

## Build a dendogram to visualize the clusters 
sampleTree = hclust(dist(datExpr0_Infect), method = "average");

## Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9);
## pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#These dataset is 88% from the original size. Zeros across all samples were removed. 
plot(sampleTree, main = "A.thaliana Infected (10 samples/threshold=17/15.5/9.39)", sub="-", xlab="Genes", ylab = "Height (log2)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

## Set a cutoff line  
abline(h = 17, col = "red");   ## 4 clust
abline(h = 15.5, col = "green"); ## 3 clust / unassigned=1
abline(h = 9.39, col="blue");     ## 24 clust / Based on the interquartile grade of the dataset (Q1-Q3)

clust = cutreeStatic(sampleTree, cutHeight=9.39, minSize = 5)
## with 0 meaning unassigned. The largest cluster is conventionally labeled 1, the next largest 2, etc.
table(clust);


