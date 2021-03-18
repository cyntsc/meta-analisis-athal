##########################################################
###
### Goal: WGCNA for data analysis of coexpression values of A.Thaliana infected by fungi
### Made by: Cynthia Soto
### Date: October 14th, 2020
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by 8 CONTROL samples of A.Thaliana H E A L T H Y  and
###                       12 samples of A. Thaliana I N F E C T  E D for ACOMYCETE FUNGY.
### 2) Row RNASeq data comes from SRA-NCBI
### 3) Data were preprocessed, alignment with STAR Tool and quantify with HTSeq-count Tool
###                         a genomic reference approach was used. 
### 4) Data were tested with several quality tests before to include in the expression matrix 
###    FastQC / HTSeq-qa / Stats: Central-Measurements & descriptive stat
### 5) Zeros at 100% across all samples were PREVIOUSLY removed  
###
### Last update: 3 December, 2020
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
athalData <- read.csv(here("data", "all_log2_dropped20ceros.csv"), header=TRUE, sep='\t')
dim(athalData);
names(athalData);
athalData[1:3]
## Some useless columns and 2 samples were removed from the dataset due their low quality 
##athalData = subset(athalData, select = -c(X,zero.counter,zero.counter)); # 24326 genes (rows) x 19 samples (cols)
athalData = subset(athalData, select = -c(X,zero.counter,SRR6283147,SRR6283148,zero.counter)); # 24326 genes (rows) x 19 samples (cols)
names(athalData);
## remove the auxiliary data and transpose the expression data for further analysis.
datExpr0 = as.data.frame(athalData[, -c(1:1)]);
#datExpr0[1:20];

# This are the sample names
names(datExpr0);  ## samples
# This are the gene names
rownames(datExpr0) = athalData$Genes;
rownames(datExpr0);

#We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts.  If not, we remove the offending genes and samples from the data

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));# Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]	
}


## Build a dendogram to visualize the clusters 
sampleTree = hclust(dist(datExpr0), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9);

# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#These dataset is 88% from the original size. Zeros across all samples were removed. 
plot(sampleTree, main = "A.thaliana Healthy and Infected (18 samples/threshold=17/15.5/9.39).", sub="", xlab="Genes", ylab = "Height (log2)", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

## Set a cutoff line  
#abline(h = 17.5, col = "red");   ## AT5G19470
abline(h = 17, col = "red");
abline(h = 15.5, col = "green");
abline(h = 9.39, col="blue");     ## 24 clust / Based on the interquartile grade of the dataset (Q1-Q3)

clust = cutreeStatic(sampleTree, cutHeight = 9.39, minSize = 5)
## with 0 meaning unassigned. The largest cluster is conventionally labeled 1, the next largest 2, etc.
table(clust);

# clust 1 contains the genes we want to keep.
keepSamples = (clust==1)
keepSamples;
datExpr = datExpr0[keepSamples, ]
nGenes = nrow(datExpr)
nSamples = ncol(datExpr)

########################### Loading meta-data: GO function annotation ontology ##############################

## metadata will help to link objects with their source sample
## geneGO_metadata will help to link objects with GO function annotation (Araport11)

traitData <- read.csv(here("data", "geneGO_metadata.csv"), sep='\t')    ## header=TRUE,
dim(traitData);
names(traitData);

# remove columns that hold information we do not need.
#allTraits = traitData[, -c(2)];
allTraits=traitData;
allTraits = allTraits[, c(2,4:5)];
allTraits;
dim(allTraits);
names(allTraits);

# Form a data frame analogous to expression data that will hold the gene function
geneGO_ID =rownames(datExpr);
geneGO_ID;
traitRows = match(geneGO_ID, allTraits$Genes);
traitRows;

datTraits = allTraits[traitRows, -1];

rownames(datTraits) = allTraits[traitRows, 1];
rownames(datTraits);

collectGarbage();

# We  now  have  the  expression  data  in  the  variable datExpr,  and  the  corresponding  clinical  traits  in  the  variable datTraits.
#  Before we continue with network construction and module detection, we visualize how the clinical traits relate to the sample dendrogram.

# Re-cluster 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors((datTraits), signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Dataset1: A.thaliana (healthy and infected) with GO annotation heatmap", 
                    cex.lab = 1, cex.axis = 1, cex.main = 1.2, cex.dendroLabels = 0.5)


##########################################################################################################
###
###                              ASSIGNING A SMALLER CLUSTERS
###
##########################################################################################################

# clust 1 contains the genes we want to keep.
keepSamples = (clust==4)
keepSamples;
datExpr = datExpr0[keepSamples, ]
nGenes = nrow(datExpr)
nSamples = ncol(datExpr)

########################### working with the meta-data: GO function annotation ontology ##############################

# Form a data frame analogous to expression data that will hold the gene function
geneGO_ID =rownames(datExpr);
geneGO_ID;
traitRows = match(geneGO_ID, allTraits$Genes);
traitRows;

datTraits = allTraits[traitRows, -1];

rownames(datTraits) = allTraits[traitRows, 1];
rownames(datTraits);

collectGarbage();

# We  now  have  the  expression  data  in  the  variable datExpr,  and  the  corresponding  clinical  traits  in  the  variable datTraits.
#  Before we continue with network construction and module detection, we visualize how the clinical traits relate to the sample dendrogram.

# Re-cluster 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors((datTraits), signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Dataset1: A.thaliana (healthy and infected) with GO annotation heatmap / Cluster 3", 
                    cex.lab = 1, cex.axis = 1, cex.main = 1.2, cex.dendroLabels = 0.5)


