#################################################################################################
###
### Goal: Analysis of Gene Significance (GS) and Module Membership (MM) for the control dataset
###
### Made by: Cynthia Soto 
### Date: Aug 17, 2022
### Latest update: xxxxx
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
###   1) Display the correlation values within a heat map plot 
###   2) Intra-modular analysis: identifying genes with high GS and MM 
###      a) create a linear regression model using R’s lm() 
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
lnames = load(file = "Athal_Healthy_MatrixD.RData");                  
# Load network data saved in the second part.
lnames = load(file = "Athal_Healthy_Module_Identification_MatrixD.RData");   
# Load de module information
lnames = load(file = "Athal_Healthy_MergedMods_MatrixD.RData");

##############     Quantifying module–trait associations     ##############     


nGenes = ncol(athalDataHealthy);
nSamples = nrow(athalDataHealthy);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(athalDataHealthy, moduleColors)$eigengenes
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
               main = paste("Module-trait healthy plants by hpi"))


###################################################################################################
##                     Intramodular analysis: identifying genes with high GS and MM
##
##                 Relating modules to external information and identifying important genes   
##                            Gene Significance (GS) and Module Membership (MM)
##                      ANALYSIS FOCUSED IN ARABIDOPSIS HEALTHY AT 30 HPI
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

# Define variable Healthy30 containing the Healthy30 column of datTrait2
Healthy30 = as.data.frame(datTraits2$healthy30);
names(Healthy30) = "Healthy30hpi"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(athalDataHealthy, MEs, use = "p"));                  # module membership MM is the correlation of the module eigengene and the gene expression profile.
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(athalDataHealthy, Healthy30, use = "p"));               # Gene Significance GS is the absolute value of the correlation between the gene and the trait.
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Healthy30), sep="");
names(GSPvalue) = paste("p.GS.", names(Healthy30), sep="");

##############  Intramodular analysis: identifying genes with high GS and MM  ##############

####    Positive correlated
modNames
#module = "coral3"

####  Negative correlated
module ="darkolivegreen4"
  
column = match(module, modNames);
moduleGenes = moduleColors==module;

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
         xlab =paste("Module Membership in", module), 
         ylab ="Gene significance for healthy_30", col = module)

# create a linear regression model using R’s lm() function and we’ll get the summary output using the summary() function.
fit <-lm(y~x)
summary(fit)
pintra = round(summary(fit)$coefficients[,4][[2]], digits=5) #p-value  
abline(fit, col='blue')
legend('topleft', col = "blue", lty = 1, box.lty = 1 ,legend = 'regression line', cex= 0.8) 
#mtext(paste('correlation = ', r, ",", "P-value = ", pintra))
mtext(paste("P-value = ", pintra))
title(paste("Module membership and Gene significance\n"))

# identifying top15 genes with high GS and MM
intra_modular_analysis = data.frame(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(athalDataHealthy)[moduleColors==module] #only the x module
#View(intra_modular_analysis)
top15_genes = intra_modular_analysis[1:15,]

# For coral3:
file_name = paste("../results-data/All_Top15_GS_GO_TAIR_Annotations/Healthy_30", module, "GS_MM_genes.csv", sep = "_");
# For darkolivegreen4:
file_name = paste("../results-data/All_Top15_GS_GO_TAIR_Annotations/Healthy_30", module, "GS_MM_genes_low.csv", sep = "_");
write.csv(top15_genes, file = file_name)
  
