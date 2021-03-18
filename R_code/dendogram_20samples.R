##########################################################
###
### Goal: WGCNA for data analysis of coexpression values of A.Thaliana
### Made by: Cynthia Soto
### Date: December 4th, 2020
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by 8 CONTROL samples of A.Thaliana H E A L T H Y  and0
###                       12 samples of A. Thaliana I N F E C T  E D for ACOMYCETE FUNGY.
### 2) Row RNASeq data comes from SRA-NCBI
### 3) Data were preprocessed, alignment with STAR Tool and quantify with HTSeq-count Tool
###                         a genomic reference approach was used. 
### 4) Data were tested with several quality tests before to include in the expression matrix 
###    FastQC / HTSeq-qa / Stats: Central-Measurements & descriptive stat
### 5) Zeros at 100% across all samples were PREVIOUSLY removed  
###
### Last update: 
##########################################################

rm(list = ls());

## This library contains all what is nedeed to use herarquical clustering
library(WGCNA);
## The here package allows you to set the top level of your project folder as “here” and to specify where things live relative to that location
library(here);   
here();
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

## Initial variables 
options(stringsAsFactors = FALSE);

## Load data
#here("files", "data", "all_log2_dropped20ceros.csv");

athalDataAll <- read.csv(here("data", "all_log2_tidy.csv"), header=TRUE, sep='\t')
# Take a quick look at what is in the data set:
dim(athalDataAll);
names(athalDataAll);

# remove samples damaged (ss30)
athalDataAll = subset(athalDataAll, select = -c(Ss30,Ss30.1));
names(athalDataAll);
# remove the auxiliary data and transpose the matrix.
datExpr0 = as.data.frame(t(athalDataAll[, -c(1:1)]));
datExpr0;

# This are the gene ids
names(datExpr0) = athalDataAll$Genes;
names(datExpr0); ## explore data
# This are the sample ids  
rownames(datExpr0) = names(athalDataAll)[-c(1:1)];
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

sampleTree = hclust(dist(datExpr0), method = "average");
sampleTree$merge;
sampleTree$height;
sampleTree$labels;
sampleTree$dist.method;
#sampleTree$height > 200;

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(width, height)
sizeGrWindow(12,9);

# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 1.2);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "A.thaliana treatments and controls (18samples).", sub="hClust/Average N=19559", xlab="Samples", ylab="Method: average", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
#abline(h = 250, col = "red");


##solo se retienen las ramas cuyo tamaño es al menos minSize.
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 50)
## with 0 meaning unassigned. The largest cluster is conventionally labeled 1, the next largest 2, etc.
table(clust);

# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nGenes;
nSamples = nrow(datExpr)
nSamples;

# The variable datExpr now contains the expression data ready for network analysis.

####################### Loading trait data ############################

traitData = read.csv("SampleTraits.csv", sep = '\t');
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData;
allTraits;
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
athal80Samples = rownames(datExpr);
athal80Samples;
traitRows = match(athal80Samples, allTraits$Sample);
traitRows;
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();

# We  now  have  the  expression  data  in  the  variable datExpr,  and  the  corresponding  clinical  traits  in  the  variable datTraits.
#  Before we continue with network construction and module detection, we visualize how the clinical traits relate to the sample dendrogram.

# Re-cluster 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits),main = "A.thaliana healthy at 70% Dendrogram and Sample ID heatmap")

# In the plot, shown in Fig. 2, white means a low value, red a high value, and grey a missing entry.The last step is to save the relevant expression and trait data for use in the next steps of the tutorial.

save(datExpr, datTraits, file = "athal1-80p-dataInput.RData")
