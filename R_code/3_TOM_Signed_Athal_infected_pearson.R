##########################################################
###
### Goal: WGCNA Std / Signed Ntw with dynamic cut off (pearson) 
###
### Method: TOM calculation: adjacency   ( UNSIGNED NTW / Pearson corr)
###
###
### Made by: Cynthia Soto 
### Date: March 01, 2021
### Latest update: xxx
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by                             ** 17     I N F E C T E D     ** 
###                                                            SAMPLES OF A.THALIANA       
###    (Raw RNA-Seq data comes from the SRA-NCBI, were aligned with STAR and were quantified with HTSeq)
### 2) All data are normalized in Log2 scale, were cleaned and filtered in different cut off thresholds
###    a) First dataset: zeros removed across samples.
###    b) Second dataset: zeros filtered dataset with cutoff in the cuartil Q1.
###    c) Third dataset: zeros filtered dataset with cutoff in the cuartil Q2.
###    d) Fourth dataset: zeros filtered dataset with cutoff in cuartil Q3.
###    ** Resulted in different sizes. 

### Some code references: 
### https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
### http://pklab.med.harvard.edu/scw2014/WGCNA.html 
### https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA


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
here();
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

## Allow multi-treads and setting this in the bash environ
allowWGCNAThreads();
ALLOW_WGCNA_THREADS=30;
## Initial variables 
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

## Load data without cut off sets

## Data just filtered of zeros ########################################################
here("files", "data", "all_infected_Log2_17ceros_drop.csv");
athalData3 <-  read.csv(here("data", "all_infected_Log2_17ceros_drop.csv"), header=TRUE,  sep='\t')

dim(athalData3);
names(athalData3);
athalData3 = subset(athalData3, select = -c(zero.counter));   # cut columns not desired

## rearrange columns accordint to traits meta-data (must match)
col_order <- c("Genes", "Bc12",   "Bc12.1", "Bc18",   "Bc18.1", "Bc24",
               "Bc24.1", "Ch22",   "Ch22.1", "Ch22.2", "Ch22.3",
               "Ch40",   "Ch40.1", "Ch40.2", "Ch40.3", "Ss30",
               "Ss30.1", "Ss30.2")
athalData3 <- athalData3[, col_order]
## You see that genes are listed in a column named "Genes" and samples are in columns
athalData3[1:5,1:10]

## Manipulate athalData3 so it matches the format WGCNA needs
row.names(athalData3) = athalData3$Genes   ## Pull the gene names 
#gene.names=rownames(athalData3)

## Transpose data to set as WGCNA needs
athalData3$Genes = NULL             ## Set genes column as index row
athalData3 = as.data.frame(t(athalData3)) # transpose dataset, now samples are rows and genes are columns
dim(athalData3) 
athalData3[1:5,1:10]
rownames(athalData3)   ## these are the samples

## Run this to check if there are gene outliers
gsg = goodSamplesGenes(athalData3, verbose = 3)
gsg$allOK

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
#if (!gsg$allOK)
#   {if (sum(!gsg$goodGenes)>0)
#       printFlush(paste("Removing genes:", paste(names(athalData3)[!gsg$goodGenes], collapse= ", ")));
#       if (sum(!gsg$goodSamples)>0)
#           printFlush(paste("Removing samples:", paste(rownames(athalData3)[!gsg$goodSamples], collapse=", ")))
#       athalData3= athalData3[gsg$goodSamples, gsg$goodGenes]
#       }

## Create an object called "datTraits" that contains the trait data
datTraits <- read.csv(here("data", "Athal_Traits_Feb2021.tsv"), header=TRUE,  sep='\t')
datTraits[1:5,]
datTraits = subset(datTraits, select = -c(accesion));
datTraits[1:5,]
## form a data frame analogous to expression data that will hold the infection traits.
rownames(datTraits) = datTraits$ID
datTraits$ID = NULL                 ## Set ID sample column as index row
table(rownames(datTraits)==rownames(athalData3)) #should return TRUE if datasets align correctly, otherwise your names are out of order

# have finished uploading and formatting expression and trait data
# Expression data is in datExpr3, corresponding traits are datTraits

save(athalData3, datTraits, file="SamplesAndTraits.RData")
load("SamplesAndTraits.RData")

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

##Plot the results
sizeGrWindow(5, 5)
par(mfrow = c(1,3));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, Signed R^2 (pearson)",
     type="n", main = paste("Scale independence. A.thal infected (Q3"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.817,col="red")
#abline(h=-0.90,col="blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=108,col="blue");

# Median connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,6],xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",main = paste("Median connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
abline(h=70.5 ,col="green");

softPower = 28;


## No hay umbral de corte optimo --> abortar !!!!


###########################################################################
#
#   Generating adjacency and TOM similarity matrices based on the selected softpower
#
###########################################################################

#calclute the adjacency matrix
adj= adjacency(athalData3,type = "signed",    #specify network type
               power = softPower);
head(adjacency)
dim(adj)
#adj[1:100]

## Construct Networks- USE A SUPERCOMPUTER -----------------------------
## translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
## this action minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(athalData3,
                          TOMType = "signed", 
                          power = softPower);
dissTOM=1-TOM

###########################################################################
#
# Generate Modules --------------------------------------------------------
#
###########################################################################

## Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM),
                      method="average");                   #"ward","single","complete","average","mcquitty","median"or"centroid"

## Plot the results
sizeGrWindow(9, 12)
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="Gene clusters", 
     main="Clustering tree based on the TOM dissimilarity measure of Athal infected (Q3)", 
     sub="Signed-Ntw (average)", cex=0.3)                  # , labels= FALSE, hang=0.04

###################################################################################
##
##    Module identification using dynamic tree cut 
##
###################################################################################

## Set the minimum module size
minModuleSize = 20;
#diff(geneTree$height)

## This chunck is the function for pruning of HClust with the method tree.
dynamicMods = cutreeDynamic(dendro = geneTree,
                            method="tree",
                            deepSplit=2,
                            #cutHeight = 0.99,
                            pamRespectsDendro= FALSE, 
                            minClusterSize = minModuleSize);

## This chunck is the function for pruning of HClust with the method disstM
dynamicMods = cutreeDynamic(dendro= geneTree, 
                            distM = dissTOM,        ## Method "hybrid". The distance matrix used as input to hclust
                            deepSplit=2,            ## For method "hybrid",  range 0 to 4. If TRUE, method will favor sensitivity and produce more smaller clusters. When FALSE, there will be fewer bigger clusters.
                            pamRespectsDendro= FALSE,     ## La etapa PAM respetará el dendrograma en el sentido de que lospequeños clusters solo se asignarán a los clusters que pertenezcan a la misma rama 
                            minClusterSize = minModuleSize)

## when cutHeight not given, for method=="tree" it defaults to 0.99, for method=="hybrid" it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram

## gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
length(dynamicMods)
sort(table(dynamicMods), decreasing = TRUE)

## Calculates module eigengenes (1st principal component) of modules in a given single dataset.
dynamicColors = labels2colors(dynamicMods)      ## colors code is assigned to plot
sort(table(dynamicColors), decreasing = TRUE)   ## Get the number of genes by color

MEList= moduleEigengenes(athalData3, 
                         colors= dynamicColors,
                         excludeGrey = TRUE,
                         softPower = softPower)
MEs= MEList$eigengenes                         ## Here are the ME tables by color module
MEDiss= 1-cor(MEs)                             ## Here is the disimmility corr of the MEs
METree= flashClust(as.dist(MEDiss), 
                   method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

## plots tree showing how the eigengenes cluster together
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="clusterwithoutmodulecolors.pdf")
plot(METree, main= "Clustering of module eigengenes for A thaliana infected (Q3)", 
     xlab= "", 
     sub= "")

## set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.0
merge = mergeCloseModules(athalData3, 
                          dynamicColors, 
                          cutHeight= MEDissThres, 
                          verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

length(mergedMEs)   ## number if merged MEs
#mergedMEs$MEviolet

#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()
#**********************************************************************************

## plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="cluster.pdf")
plotDendroAndColors(geneTree, main= "Clustering of MEs merged on 0.2 for A thaliana infected (Q3)", 
                    cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, 
                    hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")


## Write MEs in the path location
length(dynamicColors)
module_colors = setdiff(unique(dynamicColors), "grey")
length(module_colors)
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

## Write MEs merged in the path location
#length(mergedColors)
#module_colors = setdiff(unique(mergedColors), "grey")
#length(module_colors)
##module_colors
#for (color in module_colors){
#  module=SubGeneNames[which(mergedColors==color)]
#  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
#}

## Correlate traits --------------------------------------------------------

## Define number of genes and samples
nGenes = ncol(athalData3)
nSamples = nrow(athalData3)
## Recalculate MEs with color labels
MEs0 = moduleEigengenes(athalData3, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, 
                     datTraits, 
                     use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 
                                     nSamples)

## Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
#par(mar= c(6, 8.5, 3, 3))
sizeGrWindow(9, 12)
par(mar= c(6, 8.5, 3, 3))

#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.4,
               zlim= c(-1,1),
               main= paste("Module-trait relationships for A thaliana infected (Q3)"))
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

## Calculation of the topological overlap matrix

###################################################################################################

###################################################################################################
