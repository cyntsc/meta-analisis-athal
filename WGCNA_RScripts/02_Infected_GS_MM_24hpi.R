#################################################################################################
###
### Goal: Analysis of Gene Significance (GS) and Module Membership (MM) for Bc24hpi Y Ch22hpi
###
### Made by: Cynthia Soto 
### Date: Feb 18, 2022
### Latest update: May 12, 2022
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
###   1) Display the correlation values within a heatmap plot 
###   2) Intramodular analysis: identifying genes with high GS and MM 
###      a) create a linear regression model using R’s lm() 
###      b) high intramodular connectivity ~ high kwithin
###   3) Merge the statistical information of a module with gene annotation and output the summary to a file
#################################################################################################


# BASE CONFIGURATION ############################################################################

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Clear object lists
rm(list = ls());

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
library(here);    ## allows you to set the top level of your project folder as “here” and to specify where things live relative to that location   
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

# Load the expression data and trait data saved in the first part
lnames = load(file = "Athal_Infected_MatrixE.RData");                  
# Load network data saved in the second part.
lnames = load(file = "Athal_Infected_Module_Identification_MatrixE.RData");   
# Load de module information
lnames = load(file = "Athal_Infected_MergedMods_MatrixE.RData");

##############     Quantifying module–trait associations     ##############     


nGenes = ncol(athalData3);
nSamples = nrow(athalData3);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(athalData3, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

## suitable representation in a heatmap
#sizeGrWindow(10,6)
#par(mar = c(6, 8.5, 3, 3));

sizeGrWindow(9, 12)
par(mar= c(3.5, 10, 2, 1))    # margen inferior - y(izq) + margen superior - y(der)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),   # greenWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait by fungi and hpi"))


###################################################################################################
##                     Intramodular analysis: identifying genes with high GS and MM
##
##                 Relating modules to external information and identifying importantgenes   
##                            Gene Significance (GS) and ModuleMembership (MM)
##                      ANALYSIS FOCUSED IN ARABIDOPSIS WITH BOTRYTIS CINEREA A LAS 24 HPI
##
###################################################################################################

# bit refreser of the fit summary report:
# Call: This is an R feature that shows what function and parameters were used to create the model.
# Residuals: Difference between what the model predicted and the actual value of y.  
# You can calculate the Residuals section like so: summary(y-fit$fitted.values)
# Coefficients: These are the weights that minimize the sum of the square of the errors. 
# Residual std error: std.dev is the square root of variance.  Standard Error func in R is very similar. 
# Multiple R-Squared: Also called the coefficient of determination, this is an oft-cited measurement of how well your model fits to the data.
#                     it’s a quick and pre-computed check for your model.
# Adjusted R-Squared normalizes Multiple R-Squared by taking into account how many samples you have and how many variables you’re using.
# F-Statistics: this is the second “test” that the summary function produces for lm models.  The F-Statistic is a “global” test that checks if at least one of your coefficients are nonzero.
# More details: https://www.learnbymarketing.com/tutorials/explaining-the-lm-summary-in-r/ 

# Define variable Bc_24 containing the Bc_24 column of datTrait2
Bc_24 = as.data.frame(datTraits2$B_24hpi);
names(Bc_24) = "Bc_24hpi"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"));                  # module membership MM is the correlation of the module eigengene and the gene expression profile.
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(athalData3, Bc_24, use = "p"));               # Gene Significance GS is the absolute value of the correlation between the gene and the trait.
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Bc_24), sep="");
names(GSPvalue) = paste("p.GS.", names(Bc_24), sep="");

##############  Intramodular analysis: identifying genes with high GS and MM  ##############

####    Positive correlated
module = "chocolate"
module = "lavender"
module = "mistyrose3"
module = "firebrick2"

####   Negative correlated
module = "indianred"
module = "coral4"
module = "green"

####   Consensus module
module ="darkmagenta"

column = match(module, modNames);
moduleGenes = moduleColors==module;

# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for body Bc_24",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module);

#Intramodular analysis: identifying genes with high GS and MM
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
x <- abs(geneModuleMembership[moduleGenes, column])
y <- abs(geneTraitSignificance[moduleGenes, 1])
limit <- range(c(x,y)) 
r <- round(cor(x,y),2)
plot(    x, 
         y, 
         #ylim =range(min(y),1), xlim =range(min(x),1), 
         ylim =range(0,1), xlim =range(0,1), 
         xlab =paste("Module Membership in", module, "module"), 
         ylab ="Gene significance for Bc24", col = module)

# create a linear regression model using R’s lm() function and we’ll get the summary output using the summary() function.
fit <-lm(y~x)
summary(fit)
pintra = round(summary(fit)$coefficients[,4][[2]], digits=5) #p-value  
abline(fit, col='blue')
legend('topleft', col = "blue", lty = 1, box.lty = 1 ,legend = 'regression line', cex= 0.8) 
mtext(paste('correlation = ', r, ",", "P-value = ", pintra))
title(paste("Module membership vs. gene significance\n"))
# text(y~x, 
#      labels = rownames(athalData3)[moduleColors==module], 
#      data = athalData3, cex=0.8, font=4, pos = 2) #add label


# identifying genes with high GS and MM
intra_modular_analysis = data.frame(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(athalData3)[moduleColors==module] #only the x module
#View(intra_modular_analysis)
intra_modular_analysis[1:15,]

# high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of the each driver gene in the x module to all other genes in the x module)
connectivity = intramodularConnectivity(adjacency, moduleColors)
connectivity = connectivity[colnames(athalData3)[moduleColors==module],] #only the x module
order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
connectivity = connectivity[order.kWithin,] #order rows following kWithin
dim(connectivity)
#View(connectivity)
top_connectivity = connectivity[1:15,] #top 5 genes that have a high connectivity to other genes in the x module
#View(top_connectivity)
top_connectivity = connectivity[380:386,] #top 5 genes that have a high connectivity to other genes in the x module


###################################################################################################
##                     Intramodular analysis: identifying genes with high GS and MM
##
##                 Relating modules to external information and identifying importantgenes   
##                            Gene Significance (GS) and ModuleMembership (MM)
##                      ANALYSIS FOCUSED IN ARABIDOPSIS WITH COLLETOTRICHUM A LAS 22 HPI
##
###################################################################################################

# Define variable Ch_22 containing the Ch22 column of datTrait2
Ch_22 = as.data.frame(datTraits2$Ch_22hpi);
names(Ch_22) = "Ch_22hpi"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(athalData3, Ch_22, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Ch_22), sep="");
names(GSPvalue) = paste("p.GS.", names(Ch_22), sep="");

##############  Intramodular analysis: identifying genes with high GS and MM  ##############

####    High correlated
module = "chocolate2"
module = "mediumpurple1"
module = "lightsteelblue"
module = "tan3"
module = "brown2"

####    Down correlated
module = "coral4"

####   aNTAGONIC
module = "indianred"
module = "green"

####   Consens module
module ="darkmagenta"

column = match(module, modNames);
moduleGenes = moduleColors==module;

#Intramodular analysis: identifying genes with high GS and MM
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
x <- abs(geneModuleMembership[moduleGenes, column])
y <- abs(geneTraitSignificance[moduleGenes, 1])

#limit <- range(c(x,y)) 
r <- round(cor(x,y),2)
plot(    x, 
         y, 
         ylim =range(min(y),1), xlim =range(min(x), 1), 
         #ylim =range(0,1), xlim =range(0.5,1), 
         xlab =paste("Module Membership in", module, "module"), 
         ylab ="Gene significance for Ch22", col = module)

# create a linear regression model using R’s lm() function and we’ll get the summary output using the summary() function.
fit <-lm(y~x)
summary(fit)
# bit refreser:https://www.learnbymarketing.com/tutorials/explaining-the-lm-summary-in-r/ 

pintra = round(summary(fit)$coefficients[,4][[2]], digits=5) #p-value  
abline(fit, col='blue')
legend('topleft', col = "blue", lty = 1, box.lty = 1 ,legend = 'regression line', cex= 0.8) 
mtext(paste('correlation = ', r, ",", "P-value = ", pintra))
title(paste("Module membership vs. gene significance\n"))
# text(y~x, 
#      labels = rownames(athalData3)[moduleColors==module], 
#      data = athalData3, cex=0.8, font=4, pos = 2) #add label


# identifying genes with high GS and MM
intra_modular_analysis = data.frame(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(athalData3)[moduleColors==module] #only the x module
#View(intra_modular_analysis)
intra_modular_analysis[1:15,]

# high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of the each driver gene in the x module to all other genes in the x module)
connectivity = intramodularConnectivity(adjacency, moduleColors)
connectivity = connectivity[colnames(athalData3)[moduleColors==module],] #only the x module
order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
connectivity = connectivity[order.kWithin,] #order rows following kWithin
dim(connectivity)
#View(connectivity)
top_connectivity = connectivity[1:15,] #top 5 genes that have a high connectivity to other genes in the x module
#View(top_connectivity)
top_connectivity = connectivity[380:386,] #top 5 genes that have a high connectivity to other genes in the x module


############  Summary output of network analysis results   #############
##
## merge this statistical information with gene annotation and write out a
## file that summarizes the most important results and can be inspected in standard spreadsheet
##
########################################################################

# Genes to annotate included in the module
names(athalData3)                                  # these are all genes included in the analysis
names(athalData3)[moduleColors=="darkmagenta"]     # these are all genes included in the darkmagenta module

# Arabidopsis annotation file avaulable in Bioconductor 
# http://bioconductor.org/packages/2.4/Arabidopsis_thaliana.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("BSgenome.Athaliana.TAIR.04232008")

# lnames = load("Athal_Infected_AnnotFile.RData")
data_annot <- read.csv(here("data","geneGO_metadata.csv"), header = TRUE, sep = '\t');
names(data_annot)
data_annot[1:10,1:5]
annot <- data_annot[,2:5]      #subset keeping certain columns
names(annot)
annot[1:10,1:4]

# Genes del modulo que deseamos anotar 
genes = names(athalData3)
length(genes)
genes2annot = match(genes, annot$Genes)
genes2annot[1:10]
# The following is the number or genes without annotation:
sum(is.na(genes2annot))
# Should return 0.
#save(genes, annot, genes2annot, file="Athal_Infected_AnnotFile.RData")

## The modules will be ordered by their significance for <Bc_24> or <Ch_22>, with the most significant ones to the left.

# Create the starting data frame (table)
geneInfo0 = data.frame(GenesID = genes,
                       GO_term = annot$GO_term[genes2annot],
                       GO_ID = annot$GO_id[genes2annot],
                       TAIR_ID = annot$TAIR_id[genes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

length(geneInfo0)
geneInfo0[1:10,2:7]


# Order modules by their significance (GS) for Ch_22
modOrder = order(-abs(cor(MEs, Bc_24, use = "p")));

# Order modules by their significance (GS) for Ch_22
# modOrder = order(-abs(cor(MEs, Ch_22, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Bc_24: Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Bc_24hpi));
# Ch_22: Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Ch_22hpi));

geneInfo = geneInfo0[geneOrder, ]
dim(geneInfo)
length(geneInfo)
geneInfo[1:10,2:12]

write.csv(geneInfo, file = "results/geneInfo_GSignificance_Bc24hpi.csv")
# write.csv(geneInfo, file = "results/geneInfo_GSignificance_Ch22hpi.csv")
# view(geneInfo)
