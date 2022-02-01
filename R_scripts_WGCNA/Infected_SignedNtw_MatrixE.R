##########################################################
###
### Goal: WGCNA Std / Signed Ntw with dynamic cut off (pearson) 
###
### Method: TOM calculation: adjacency   ( SIGNED NTW / Pearson corr)
###
###
### Made by: Cynthia Soto 
### Date: March 23, 2021
### Latest update: May 21, 2021
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) RNA-Seq libraries are                             ** 14  I N F E C T E D     **  w/o ATYPICAL SAMPLE DISTRIBUTIONS 
###                                                       SAMPLES OF A.THALIANA       
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
library(tidyverse);
library(dynamicTreeCut);    ## contains all what is nedeed to use herarquical clustering
library(flashClust);
library(WGCNA);
library(lattice, survival);     ## Let you calculate percentiles and more ... perhaps is required "Formula" pkg
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

#athalData3 <-  read.csv(here("data", "all_infected_Log2TPM_No_Ss30.csv"), header=TRUE,  sep='\t')
athalData3 <-  read.csv(here("data", "matrix_E_infected.csv"), header=TRUE,  sep='\t')
dim(athalData3);   # 20265  x  15

## rearrange columns accordint to traits meta-data (must match)
col_order <- c("Genes", "Bc12",   "Bc12.1", "Bc18",   "Bc18.1", "Bc24",
               "Bc24.1", "Ch22",   "Ch22.1", "Ch22.2", "Ch22.3",
               "Ch40",   "Ch40.1", "Ch40.2", "Ch40.3")
athalData3 <- athalData3[, col_order]
## You see that genes are listed in a column named "Genes" and samples are in columns
athalData3[1:5,]

## Manipulate athalData3 so it matches the format WGCNA needs
row.names(athalData3) = athalData3$Genes   ## Pull the gene names 
gene.names=rownames(athalData3)
## Pull the names according with the size grp
SubGeneNames=gene.names;

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

# TRAITMENTS BY SAMPLE (fungi-type-hpi interaction)
datTraits <- read.csv(here("data", "Traits_Infected_sample.csv"), header=TRUE,  sep='\t')   
datTraits = subset(datTraits, select = -c(accesion, Ss30, Ss30.1, Ss30.2, B_cinerea, C_higginsianum, S_sclerotium));
datTraits[1:5,]
datTraits <- datTraits[1:(nrow(datTraits)-3),]
dim(datTraits)
## form a data frame analogous to expression data that will hold the infection traits.
rownames(datTraits) = datTraits$ID
datTraits$ID = NULL                 ## Set ID sample column as index row
table(rownames(datTraits)==rownames(athalData3)) #should return TRUE if datasets align correctly, otherwise your names are out of order

# TRAITMENTS BY HPI
datTraits2 <- read.csv(here("data", "Traits_Infected_hpi.csv"), header=TRUE,  sep=',')   
datTraits2 = subset(datTraits2, select = -c(accesion, Ss_30hpi, B_cinerea, C_higginsianum, S_sclerotium));
datTraits2[1:5,]
datTraits2 <- datTraits2[1:(nrow(datTraits2)-3),]
dim(datTraits2)
## form a data frame analogous to expression data that will hold the infection traits.
rownames(datTraits2) = datTraits2$ID
datTraits2$ID = NULL                 ## Set ID sample column as index row
table(rownames(datTraits2)==rownames(athalData3)) #should return TRUE if datasets align correctly, otherwise your names are out of order

save(athalData3, datTraits, datTraits2, file="Athal_Infected_MatrixE.RData")
#load("Athal_Infected_MatrixE.RData")

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
#write.table(sft, file = "results/SFT_corr_Athal_Infected_MatrixE.txt", sep = ",", quote = FALSE, row.names = T)

##Plot the results
sizeGrWindow(5, 5)
par(mfrow = c(1,3));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, Signed R^2 / Log2(TPM)",
     type="n", main = paste("Scale independence for Athal.infected \n (Matrix E)"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red");  

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=271,col="blue");

# Median connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,6],xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",main = paste("Median connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
abline(h=207 ,col="green");

softPower = 27; # all_infected NO Ss30 samples

###########################################################################
#
#   Generating adjacency and TOM diss/similarity matrices based on the selected softpower
#
###########################################################################

#calclute the adjacency matrix
adjacency= adjacency(athalData3,type = "signed",    #specify network type
               power = softPower);
head(adjacency)
dim(adjacency)
adj[1:20]

## translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity
## this action minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(adjacency,                         #athalData3
                          TOMType = "signed", 
                          power = softPower);
dissTOM=1-TOM

save(dissTOM, file="Athal_Infected_DissTOM_MatrixE.RData")
load("Athal_Infected_DissTOM_MatrixE.RData")

###########################################################################
##
##    Generate Modules 
##
###########################################################################

## Generate a clustered gene tree
#geneTree = flashClust(as.dist(dissTOM),method="average")      #"ward","single","complete","average","mcquitty","median"or"centroid"
geneTree = hclust(as.dist(dissTOM), method="average")

## Plot the results
sizeGrWindow(9, 12)
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="Gene clusters", main="Gene clustering on TOM−based dissimilarity for A.thaliana infected (Matrix E)", cex=0.3)   # sub="Signed-Ntw (average)" , labels= FALSE, hang=0.04

###################################################################################
##
##    Module identification using dynamic tree cut 
##
###################################################################################

## Set the minimum module size
minModuleSize = 20;
#diff(geneTree$height)

## Method1: This chunck is the function for pruning of HClust with the method tree.
# dynamicMods = cutreeDynamic(dendro = geneTree,
#                            method="tree",
#                            deepSplit=TRUE,
#                            #cutHeight = 0.99,
#                            pamRespectsDendro= FALSE, 
#                            minClusterSize = minModuleSize);

## Method2: This chunck is the function for pruning of HClust with the method disstM
## We chose dissTOM for this analysis because all genes are asigned to some module *** ###
dynamicMods = cutreeDynamic(dendro= geneTree, 
                            distM = dissTOM,        ## Method "hybrid". The distance matrix used as input to hclust
                            deepSplit=2,            ## For method "hybrid",  range 0 to 4. If TRUE, method will favor sensitivity and produce more smaller clusters. When FALSE, there will be fewer bigger clusters.
                            pamRespectsDendro= FALSE,     ## La etapa PAM respetará el dendrograma en el sentido de que los pequeños clusters solo se asignarán a los clusters que pertenezcan a la misma rama 
                            minClusterSize = minModuleSize)

## when cutHeight not given, for method=="tree" it defaults to 0.99, for method=="hybrid" it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram

## gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)
# write.table(table(dynamicMods), file = "results/dynamicMods_Athal_Infected_MatrixE.txt", sep = ",", quote = FALSE, row.names = F)

## # Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)      ## colors code is assigned to plot
sort(table(dynamicColors), decreasing = TRUE)   ## Get the number of genes by color
# write.table(sort(table(dynamicColors), decreasing = TRUE), file = "results/dynamicColors_Athal_Infected_MatrixE.txt", sep = ",", quote = FALSE, row.names = F)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for Athal infected \n Matrix E")

save(dynamicMods, dynamicColors, file="Athal_Infected_DynamicMods_MatrixE.RData")
load("Athal_Infected_DynamicMods_MatrixE.RData")

##################################################################################################
#
##                      Calculation of the topological overlap matrix
#
##################################################################################################
## Raise the dissimilarity matrix to the power of 27 to bring out the module structure
sort(table(dynamicColors), decreasing = TRUE)
restGenes= (dynamicColors != "grey")                   #restGenes= (dynamicColors == "mediumpurple2")
length(restGenes)
##set the DIAGONAL of the dissimilarity to NA 
diag(dissTOM) = NA;
sizeGrWindow(7,7)
TOMplot(dissTOM, 
        geneTree, 
        as.character(dynamicColors[dynamicColors]),
        main = "TOM for A.thaliana infected / Matrix E");

##################################################################################################
#
#                          Merging of modules whose expression profiles are very similar
#
##################################################################################################

# Calculate eigengenes
MEList= moduleEigengenes(athalData3, colors= dynamicColors,
                         excludeGrey = TRUE,
                         softPower = softPower)
MEs = MEList$eigengenes                         ## Here are the ME tables by color module
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)          
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss),method= "average")
save(MEList, MEs, MEDiss, METree, file= "Athal_Infected_Module_Identification_MatrixE.RData")
#load("Athal_Infected_Module_Identification_MatrixE.RData")

# Clustering of module eigengenes
sizeGrWindow(6,14)
## plots tree showing how the eigengenes cluster together
#plot(METree, main = "Clustering of module eigengenes for A.thaliana infected", xlab = "", sub = "", labels = FALSE)
plot(METree, main= "Clustering of module eigengenes for A.thaliana infected (Matrix E)", 
     xlab= "Module Eigengenes", 
     ylab = "", sub= "",
     cex.main = 1, cex.lab = 1, cex.axis = 1)
     
##########################################################################################
##
## set a threhold for merging modules. Set the threshold in the MEDissThres=0.0 variable
##
##########################################################################################

# We choose a height cut of 0.10, corresponding to correlation of 0.90, to merge

MEDissThres = 0.10
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(athalData3, dynamicColors, 
                          cutHeight= MEDissThres, 
                          verbose =3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs
length(mergedMEs)   ## number if merged MEs
sort(table(mergedColors), decreasing = TRUE)
# write.table(sort(table(mergedColors), decreasing = TRUE), file = "results/mergedColors_Athal_Infected_MatrixE.txt", sep = ",", quote = FALSE, row.names = F)

#mergedMEs$MEviolet
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="clusterME.pdf")

## plot dendrogram with module colors below it
plotDendroAndColors(geneTree, main= "Clustering of MEs corresponding to correlation of 0.90 merged \nA thaliana infected (Matrix E)", 
                    cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, 
                    hang=0.03, addGuide= TRUE, guideHang=0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));    # standardColors(50)
moduleLabels = match(mergedColors, colorOrder)-1     ## moduleColors
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

save(MEs, moduleLabels, moduleColors, geneTree, file= "Athal_Infected_MergedMods_MatrixE.RData")
#load("Athal_Infected_MergedMods_MatrixE.RData")

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
  write.table(module, paste("mod_merged_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

##################################################################################################
#
##Calculation of the topological overlap matrix for merged modules 
#
##################################################################################################
## Raise the dissimilarity matrix to the power of 27 to bring out the module structure
sort(table(mergedColors), decreasing = TRUE)
##set the DIAGONAL of the dissimilarity to NA 
#diag(dissTOM) = NA;
sizeGrWindow(7,7)
TOMplot(dissTOM, 
        geneTree, 
        as.character(mergedColors[mergedColors]),
        main = "TOM for A.thaliana infected Merged at 0.10 (Matrix E)");

# ##################################################################################################
# #
# ## classical multi-dimensional scaling plots
# ## for visualizing the network. Here we chose 2 scaling dimensions
# #
# ##################################################################################################
# You should know that most of the functions require a distance matrix as the main argument to work with.
# The most popular function to perform a classical scaling is cmdscale() which comes with the default distribution of R. 

#cmd1=cmdscale(as.dist(dissTOM),3);
#par(mfrow=c(1,1));
#plot(cmd1, col=as.character(mergedColors), main='MDS plot for A.thaliana infected')

#####################################################################################
##
##       Correlate traits (Quantifying module–trait associations)
##       Identify modules that are significantly associated with the measured biology traits
##       Here, We use a correlation by sample
##
#####################################################################################

## Define the number of genes and samples
nGenes = ncol(athalData3)
nSamples = nrow(athalData3)
## Recalculate MEs with color labels.  moduleColors was replaced by mergedColors (line.317)    
MEs0 = moduleEigengenes(athalData3, moduleColors)$eigengenes  # Calculates module eigengenes (1st principal component) of modules in a given single dataset
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)   #Calculates Student asymptotic p-value for given correlations.
                                                                 #Asymptotic p-values are useful for large sample sizes when the calculation of an exact p-value is too computer-intensive.
# We color code each association by the correlation value

## Print correlation and p-value heatmap between modules and traits
## textMatrix= paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep= "")
## Print just the correlation heatmap between modules and traits.
textMatrix= paste(signif(moduleTraitCor, 2))
## Print "" if too many modules
#textMatrix= ""
dim(textMatrix) = dim(moduleTraitCor)
#par(mar= c(6, 8.5, 3, 3))
sizeGrWindow(9, 12)
par(mar= c(3, 10, 2, 1))    # margen inferior - y(izq) + margen superior - y(der)

#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels = names(MEs),
               font.lab.x = 2,                          # set the font blod style
               font.lab.y = 0.3,
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50), #blueWhiteRed(50), greenBlackRed()
               textMatrix= (textMatrix),
               setStdMargins= FALSE,
               cex.text= 0.6,
               zlim= c(-1,1),
               main= paste("Module-trait relationships by sample for A thaliana infected (Matrix E)"))


#####################################################################################
##
##       Correlate traits (Quantifying module–trait associations)
##       Identify modules that are significantly associated with the measured biology traits
##       *** Here we use correlation by HPI 
##
#####################################################################################

## Define the number of genes and samples
nGenes = ncol(athalData3)
nSamples = nrow(athalData3)
## Recalculate MEs with color labels.  moduleColors was replaced by mergedColors (line.317)    
MEs0 = moduleEigengenes(athalData3, moduleColors)$eigengenes  # Calculates module eigengenes (1st principal component) of modules in a given single dataset
MEs = orderMEs(MEs0)

moduleTraitCor2 = cor(MEs, datTraits2, use= "p")
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamples)   #Calculates Student asymptotic p-value for given correlations.
#Asymptotic p-values are useful for large sample sizes when the calculation of an exact p-value is too computer-intensive.
# We color code each association by the correlation value

## Print correlation and p-value heatmap between modules and traits
textMatrix2= paste(signif(moduleTraitCor2, 2), " (", signif(moduleTraitPvalue2, 1), ")", sep= "")
## Print just the correlation heatmap between modules and traits.
#textMatrix= paste(signif(moduleTraitCor, 2), " / (", signif(moduleTraitPvalue, 1), ")", sep= "")
## Print "" if too many modules
#textMatrix= ""
dim(textMatrix2) = dim(moduleTraitCor2)
#par(mar= c(6, 8.5, 3, 3))
sizeGrWindow(9, 12)
par(mar= c(3.5, 10, 2, 1))    # margen inferior - y(izq) + margen superior - y(der)

#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor2,
               xLabels= names(datTraits2),
               yLabels = names(MEs),
               font.lab.x = 2,                          # set the font blod style
               font.lab.y = 0.3,
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50), #blueWhiteRed(50), greenBlackRed()
               textMatrix= (textMatrix2),
               setStdMargins= FALSE,
               cex.text= 0.6,
               zlim= c(-1,1),
               main= paste("Module-trait by hpi for A thaliana infected (Matrix E)"))



















###################################################################################################
##
##                 Relating modules to external information and identifying importantgenes   
##                                   Gene Significance (GS) and ModuleMembership (MM)
##                 Important Note: 
##                                  Chunck below is just an example
##                                  go to the file GeneSignificance_ModuleMembership.R to see more details 
###################################################################################################

# Define variable hpi containing the weight column of datTrait
hpi12 = as.data.frame(datTraits2$X12hpi);
names(hpi12) = "X12hpi"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(athalData3, hpi12, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(hpi12), sep="")
names(GSPvalue) = paste("p.GS.", names(hpi12), sep="")


## Intramodular analysis using the GS and MM values
module = "blanchedalmond"
column = match(module, modNames);
moduleGenes = moduleColors == module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));

lng=length(abs(geneModuleMembership[moduleGenes, column]))
title = paste("Module membership vs. gene significance for A. thaliana infected \n(", lng, " genes )");
x_lab = paste("Gene significance for", names(hpi12));

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, " module"),
                   ylab = x_lab,
                   main = title,
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module)


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


