##########################################################
###
### Goal: WGCNA Std / Signed Ntw with dynamic cut off (pearson) 
###
### Method: TOM calculation: adjacency   ( UNSIGNED NTW / Pearson corr)
###
###
### Made by: Cynthia Soto 
### Date: Feb 25, 2021
### Latest update: March 03, 2021
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by                             ** 17     I N F E C T E D     ** 
###                                                            SAMPLES OF A.THALIANA       
###
###    (Raw RNA-Seq data was extracted from the SRA-NCBI DB. Once ready to analysis datasets were aligned with STAR and quantified with HTSeq)
### 2) Data are normalized into Log2 scale.
###    a) Input dataset were merged and zeros across sameples were removed, dataframe "base" for subsequent analysis. 
###    b) Second dataset: Cutoff in the cuartil Q1.
###    c) Third dataset: Cutoff in the cuartil Q2.
###    d) Fourth dataset: Cutoff in cuartil Q3.

### Here some code references, but were used many more: 
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

## Load data: we have datasetes to test (*_Q1, *_Q2 and *_Q3)

## Data filtered on cuartil Q1 ########################################################
#here("files", "data", "all_infected_Log2_17ceros_Q1_drop.csv");
#athalData3 <- read.csv(here("data", "all_infected_Log2_17ceros_Q2_drop.csv"), header=TRUE, row.names='Genes', sep='\t')

## Data filtered on cuartil Q2 ########################################################
#here("files", "data", "all_infected_Log2_17ceros_Q2_drop.csv");
#athalData3 <- read.csv(here("data", "all_infected_Log2_17ceros_Q2_drop.csv"), header=TRUE, row.names='Genes', sep='\t')

## Data filtered on cuartil Q3 ########################################################

here("files", "data", "all_infected_Log2_17ceros_Q3_drop.csv");
athalData3 <-  read.csv(here("data", "all_infected_Log2_17ceros_Q3_drop.csv"), header=TRUE,  sep='\t')

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
datTraits <- read.csv(here("data", "Athal_Traits_complete_March2021.tsv"), header=TRUE,  sep='\t')   ## This is a trait arrangment by fungus: Athal_Traits_Feb2021.tsv
datTraits[1:17,]
datTraits = subset(datTraits, select = -c(accesion, B_cinerea, C_higginsianum, S_sclerotium));
datTraits[1:17,]
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
abline(h=74 ,col="green");

softPower = 27;

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
##
##    Generate Modules 
##
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
                            deepSplit=TRUE,
                            #cutHeight = 0.99,
                            pamRespectsDendro= FALSE, 
                            minClusterSize = minModuleSize);

## ***** dissTOM is the method we choose for this analysis because all genes are asigned to some module *** ###

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
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_Signed_Module_Identification.RData")
load("Network_Signed_Module_Identification.RData")

## plots tree showing how the eigengenes cluster together
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="clusterwithoutmodulecolors.pdf")
plot(METree, main= "Clustering of module eigengenes for A thaliana infected (Q3)", 
     xlab= "", 
     sub= "")


## plots tree showing how the eigengenes by color coder (add.cyntsc)
plotDendroAndColors(geneTree, main= "Clustering of MEs for A thaliana infected (Q3)", 
                    dynamicColors, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang=0.03, addGuide= TRUE, guideHang=0.05)

##########################################################################################
##
## set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
##
##########################################################################################

MEDissThres = 0.0
merge = mergeCloseModules(athalData3, 
                          dynamicColors, 
                          cutHeight= MEDissThres, 
                          verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

length(mergedMEs)   ## number if merged MEs
sort(table(mergedColors), decreasing = TRUE)

mergedMEs$MEviolet

#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()
#**********************************************************************************

## plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="cluster.pdf")
plotDendroAndColors(geneTree, main= "Clustering of MEs merged on 0.0 for A thaliana infected (Q3)", 
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

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_merge.RData")
load("Network_signed_merge.RData")


## Write modules no merged in the path location
length(dynamicColors)
module_colors = setdiff(unique(dynamicColors), "grey")    ## exclude the grey module reserved for not assigned genes
length(module_colors)
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

## Write MEs merged in the path location
length(mergedColors)
module_colors = setdiff(unique(mergedColors), "grey")     ## exclude the grey module reserved for not assigned genes
length(module_colors)
#module_colors
for (color in module_colors){
  module=SubGeneNames[which(mergedColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

#####################################################################################
##
##       Correlate traits 
##
#####################################################################################

## Define number of genes and samples
nGenes = ncol(athalData3)
nSamples = nrow(athalData3)
## Recalculate MEs with color labels
MEs0 = moduleEigengenes(athalData3, moduleColors)$eigengenes   #Calculates module eigengenes (1st principal component) of modules in a given single dataset
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, 
                     datTraits, 
                     use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor,           #Calculates Student asymptotic p-value for given correlations.
                                     nSamples)                 #Asymptotic p-values are useful for large sample sizes when the calculation of an exact p-value is too computer-intensive.

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
               cex.text= 0.7,
               zlim= c(-1,1),
               main= paste("Module-trait relationships for A thaliana infected (Q3)"))
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

## Calculation of the topological overlap matrix

###################################################################################################
##
##                 Relating modules to external information and identifying importantgenes   
##                                   Gene Significance (GS) and ModuleMembership (MM)
##                 Important Note: 
##                                  Chunck below is just an example
##                                  go to the file GeneSignificance_ModuleMembership.R to see more details 
###################################################################################################

# Define variable Bc24 containing the Bc24 column of datTrait

BotrytisC24 = as.data.frame(datTraits$Bc24)
names(BotrytisC24) = "BotrytisC24"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(athalData3, BotrytisC24, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(BotrytisC24), sep="")
names(GSPvalue) = paste("p.GS.", names(BotrytisC24), sep="")

## Intramodular analysis using the GS and MM values

module = "grey60"
column = match(module, modNames);
moduleGenes = moduleColors == module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));

title = paste("Module membership vs. gene significance for", module," module\n")
x_lab = paste("GS for A.thaliana with", names(BotrytisC24))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = x_lab,
                   main = title,
                   cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = module)

###################################################################################################
##
##                 Summary output of network analysis results
##                                
##                 
###################################################################################################

## The data expression are only annotated by ID names
names(athalData3)
## will return all probe IDs included in the analysis
length(names(athalData3)[moduleColors=="blue"])

## load a file to make it more interpretable
## Fields in the file:  Genes / GO_term /	GO_id /	TAIR_id	/ GO_slim
annot = read.csv(here("data", "geneGO_metadata.csv"), header=TRUE, sep='\t')
annot[1:5,]
dim(annot)
names(annot)
probes = names(athalData3)
#####probes = names(athalData3)[moduleColors=="blue"]
probes2annot = match(probes, annot$Genes)

## The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0


# Create the starting data frame
geneInfo0 = data.frame(gene_ID = probes,
                       GO_term = annot$GO_term[probes2annot],
                       GO_slim = annot$GO_slim[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
nrow(geneInfo0)

# Order modules by their significance for BotrytisC24
modOrder = order(-abs(cor(MEs, BotrytisC24, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.BotrytisC24));
geneInfo = geneInfo0[geneOrder, ];

view(geneInfo)

write.csv(geneInfo, file = "geneInfo_GS_MM.csv")


