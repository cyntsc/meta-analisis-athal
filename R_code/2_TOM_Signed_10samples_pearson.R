##########################################################
###
### Goal: WGCNA Std / Signed Ntw with dynamic cut off (pearson) 
###
### Method: TOM calculation: adjacency   ( SIGNED NTW / Pearson corr)
###
###
### Made by: Cynthia Soto 
### Date: January 11, 2021 / Last update: Jan 30, 2021
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Dataset is composed by 10 samples of A. Thaliana: ** I N F E C T E D ** 
### 2) Row RNASeq data comes from SRA-NCBI
###    Data were cleaned, aligned with STAR and quantified with HTSeq-count Tool
### 4) Data were tested with several quality tests before to include in the expression matrix 
###    FastQC / HTSeq-qa / Stats: Central-Measurements & descriptive stat
### 5) 100% ceros across all samples were removed  

### Reference: http://pklab.med.harvard.edu/scw2014/WGCNA.html 

##########################################################

#getwd();        ## setwd('../../data/run/cyntsc');
rm(list = ls());
# Clear plots
if(!is.null(dev.list())) dev.off()


#install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
#source("http://bioconductor.org/biocLite.R") 
#biocLite("impute")

library(tidyverse)
## This library contains all what is nedeed to use herarquical clustering
library(dynamicTreeCut)
library(flashClust)
library(WGCNA);
## This library needs lattice, survival and Formula. Let you calculate percentiles and more ...
library(lattice, survival, Formula);
library(Hmisc);   
## The here package allows you to set the top level of your project folder as “here” and to specify where things live relative to that location
library(here);   
here();
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

## Allow multi-treads and setting this in the bash environ
allowWGCNAThreads();
ALLOW_WGCNA_THREADS=30;
## Initial variables 
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

## Load data
here("files", "data", "all_log2_tidy.csv");

## 20 samples are included
athalData3 <- read.csv(here("data", "all_log2_tidy.csv"), header=TRUE, row.names='Genes', sep='\t')
dim(athalData3);
names(athalData3);

## get means for variables in data frame athalData3, excluding missing values 
sapply(athalData3, mean, na.rm=TRUE)
sapply(athalData3, range, na.rm=TRUE)
sapply(athalData3, sd, na.rm=TRUE)
sapply(athalData3, quantile, na.rm=TRUE)
summary(athalData3)

## Samples of interest are filtered: 2 samples are removed due does not match the expected % of read alignment.  
athalData3 = subset(athalData3, select = -c(Ss30,Ss30.1)); # 24326 genes (rows) x 19 samples (cols)
athalData3 = subset(athalData3, select =c(9:18)); 
dim(athalData3);   ## 10 samples are keept after applying these filters.
#class(athalData3);  #Check data object type: is a data.frame

# A BIT OF STATS: plot median and std.dev lines ********************************************************************************
stats_infected = describe(athalData3, digits=3, tabular=FALSE, scroll=TRUE)           #It is a generalization of SAS UNIVARIATE
stats_infected[sort(names(stats_infected))]
#stats_infected$Ch40;
#stats_infected[c('Ch40','Ch40.1','Ch22')]

samples=names(athalData3)
stats_infec_x = (apply(athalData3, 2, mean, na.rm=T))
sort(stats_infec_x)
stats_infec_m = (apply(athalData3, 2, median, na.rm=T))
sort(stats_infec_m)
stats_infec_sd=((sapply(athalData3, sd)))
sort(stats_infec_sd)
plot(x=seq(samples), y=stats_infec_x, type="l", lty=1, col="green", lwd=2, ylim=c(0,10), xlim=c(0,10),
     axes=F,bty="n",yaxs="i",xaxs="i", main="Basic stats of A.thaliana infected",     
     xlab="Samples", ylab="Log2 normalized values")
# plot dashed line
lines(x=seq(samples), y=stats_infec_m, lty=2, col="brown", lwd=2,)
lines(x=seq(samples), y=stats_infec_sd, lty=3, col="red", lwd=2,)
# add axes
axis(side=1, labels=samples, at=seq(samples))
axis(side=2, at=seq(0,10), las=0)
legend(0.2,10,c("mean","median","std.dev"), lwd=c(2,2,2), col=c("green","brown","red"), y.intersp=1)


## WGCNA requires genes be given in the columns ***************************************************************

## double check evaluation to avoid the error caused by empty correlation or NANs
## this code build a matrix with all genes with ceros
# athalData3 == 0
# rowSums(athalData3 == 0)                                        #sum the number of ceros per genes
# rowSums(athalData3 == 0) >=10
# sum(rowSums(athalData3 == 0) >=10)  
# athal_withceros=athalData3[rowSums(athalData3 == 0) >=10, ]     #delete rows (genes) that sum cero

## this code preserves a matrix with certain number of genes with ceros
rowSums(athalData3 == 0) < 1
sum(rowSums(athalData3 == 0) < 1)  
athalData3=athalData3[rowSums(athalData3 == 0) < 1, ]     #delete rows (genes) that sum cero

dim(athalData3)
sum(rowSums(athalData3 == 0) >= 1) 

## Pull the names and transpose data 
gene.names=rownames(athalData3)
#gene.names
athalData3=as.data.frame(t(athalData3))
athalData3[1:10,1:5]

## Focus in a small group of genes to test the function **************** JUMP IF DATA READY **********************
# n=5000
# ## transpose data
# athalData3=as.data.frame(t(athalData3))
# athalData3=athalData3[1:n]
# ## Pull the names according with the size grp
# SubGeneNames=gene.names[1:n]
## ***************************************************************************************************************

datExpr=athalData3
any(is.na(datExpr))

###########################################################################################################
#
#                    pickSoftThreshold function offers a analysis of scale free topology for soft-thresholding 
#
##########################################################################################################

## https://rdrr.io/cran/WGCNA/man/pickSoftThreshold.html 
## Choosing a soft-threshold to fit a scale-free topology to the network
powers = c(c(1:10), seq(from = 12, to=20, by=2));
powers = c(c(1:10), seq(from = 12, to=30, by=2));
## remember that a signed ntw preserve the natural continuity of the correlation (+ / -), contrary to whats an unsigned ntw does
sft=pickSoftThreshold(
  datExpr,
  dataIsExpr = TRUE,
  #RsquaredCut = 0.85,                        # desired minimum scale free topology fitting index R^2.
  powerVector = powers,
  corFnc = cor,                              # cor: Fast calculations of Pearson correlation
  # corFnc (use = "all.obs",                    #  c("pearson", "kendall", "spearman")
  #        method = c("pearson"),
  #        weights.x = NULL),
  corOptions = list(use = 'p'),             # Almost all lists in R internally are Generic Vectors, whereas traditional dotted pair lists (as in LISP) remain available but rarely seen by users (except as formals of functions). 
  networkType = "signed");                  # "unsigned", "signed", "signed hybrid"

#warnings()
sft

##Plot the results
sizeGrWindow(5, 5)
par(mfrow = c(1,3));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2 (pearson)",type="n", main = paste("Scale independence. 10 samples."));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.71,col="red")
#abline(h=-0.90,col="blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=235,col="blue");

# # Median connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,6],xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",main = paste("Median connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
# abline(h=1860,col="green");

softPower = 30;

###########################################################################
#
#   Generating adjacency and TOM similarity matrices based on the selected softpower
#
###########################################################################

#calclute the adjacency matrix
adj= adjacency(datExpr,type = "unsigned", 
               power = softPower);
dim(adj)
#adj[1:100]

##CSC: some validation to avoid the error caused by NANs. Be careful to check your data before to continue
# any(is.na(adj))                   #check for any NANs
# sum(is.na(adj))                   #sum the number of NANs
# adj[1:100]
# any(is.infinite(adj))
# any(is.null(adj))
# # na.omit(datExpr) helps in omitting NA

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,
                          networkType = "unsigned", 
                          TOMType = "unsigned", 
                          power = softPower);

# any(is.na(TOM))   #CSC
# TOM[1:100]
# dim(TOM)

##Pull genes names to the TOM matrix 
SubGeneNames=gene.names
colnames(TOM)=rownames(TOM)=SubGeneNames;
dissTOM=1-TOM

###########################################################################
#
#     Module detection
#
###########################################################################

#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),
                      method="average");                   #"ward","single","complete","average","mcquitty","median"or"centroid"

# Plot the results
sizeGrWindow(9, 12)
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="Gene clusters", main="HClust of the genes based on the TOM dissimilarity measure", sub="Method:average", cex=0.3)

# Set the minimum module size
minModuleSize = 20;

## Module identification using dynamic tree cut *******************************

diff(geneTree$height)

## Function for pruning of Hierarchical Clustering Dendrograms
dynamicMods = cutreeDynamic(dendro = geneTree,
                            method="tree",
                            cutHeight = 0.99,
                            minClusterSize = minModuleSize);
dynamicMods = cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            method="hybrid",
                            deepSplit = 2,
                            pamRespectsDendro = TRUE,          #the PAM stage will respect the dendrogram in the sense that objects and small clusters will only be assigned to clusters that belong to the same branch that the objects or small clusters being
                            # labelUnlabeled==TRUE,
                            minClusterSize = minModuleSize,
                           );
## when cutHeight not given, for method=="tree" it defaults to 0.99, for method=="hybrid" it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram

#**********************************************************************************

##gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
sort(table(dynamicMods), decreasing = TRUE)

##Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
sort(table(dynamicColors), decreasing = TRUE)

plotDendroAndColors(geneTree, 
                    dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

##discard the unassigned genes, and focus on the rest ******************************************************************

restGenes= (dynamicColors != "grey")                   #restGenes= (dynamicColors == "mediumpurple2")
length(restGenes)

###########################################################################
#
##Calculation of the topological overlap matrix
#
###########################################################################

diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], 
                              corType = "pearson",                             #"pearson" and "bicor"
                              networkType = "unsigned", 
                              power = softPower)

colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )                            #flashClust is the same that hclust but faster

plotDendroAndColors(hier1, 
                    dynamicColors[restGenes], "Dynamic Tree Cut", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors (cor.type=pearson)")    #sub="Discarded the unassigned genes"))

##consult the module's color codes
standardColors();                               #with NULL all (approx. 450) colors will be returned

##set the DIAGONAL of the dissimilarity to NA ******************************
diag(diss1) = NA;

##Visualize the TOM plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7,7)
TOMplot(diss1, 
        hier1, 
        as.character(dynamicColors[restGenes]));

###########################################################################
#
## Extract modules 
#
###########################################################################
length(dynamicColors)
module_colors= setdiff(unique(dynamicColors), "grey")
length(module_colors)
module_colors= unique((dynamicColors), "turquoise")

##module_colors
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
#module

module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
heatmap(t(m),
        zlim=c(0,1),
        col=gray.colors(100),
        Rowv=NA,
        Colv=NA,
        labRow=NA,
        scale="none",
        RowSideColors=dynamicColors[module.order])

#We can now look at the module gene listings and try to interpret their functions .. for instance using http://amigo.geneontology.org/rte

###########################################################################
#
#  Quantify module similarity by eigengene correlation. 
# Eigengenes: Module representatives 
#
###########################################################################

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))


