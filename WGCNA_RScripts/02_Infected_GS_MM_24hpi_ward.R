
# https://github.com/huynguyen250896/drivergene/blob/master/Data_n_SourceCode/2.Association-significantgenes-and-hubgenes.R
# FILE: 2.Association-significantgenes-and-hubgenes.R
# Replicating for module selection based on GS and MM

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

moduleTraitCor2 = cor(MEs, datTraits2, use= "p")   #"pearson", "kendall", "spearman")
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamples)   #Calculates Student asymptotic p-value for given correlations.
#Asymptotic p-values are useful for large sample sizes when the calculation of an exact p-value is too computer-intensive.
# We color code each association by the correlation value

## Print in a heatmap correlation and p-value between modules and traits
textMatrix2 = paste(signif(moduleTraitCor2, 2), " (", signif(moduleTraitPvalue2, 1), ")", sep= "")
## Print just the correlation heatmap between modules and traits.
#textMatrix2 = paste(signif(moduleTraitCor2, 2), " / (", signif(moduleTraitPvalue, 1), ")", sep= "")
#textMatrix2 = paste(signif(moduleTraitCor2, 2))
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
               main= paste("Module-trait by h.p.i"))


###################################################################################################
##                     Intramodular analysis: identifying genes with high GS and MM
##
##                 Relating modules to external information and identifying importantgenes   
##                            Gene Significance (GS) and ModuleMembership (MM)
##                      ANALYSIS FOCUSED IN THE SAMPLES for COLLETROTRICHUM A LAS 22 HPI
##
###################################################################################################

# Define variable Ch_22 containing the Ch22 column of datTrait2
Ch_22 = as.data.frame(datTraits2$Ch_22hpi);
names(Ch_22) = "Ch_22hpi"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"));         # module membership MM is the correlation of the module eigengene and the gene expression profile.
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(athalData3, Ch_22, use = "p"));      # Gene Significance GS is the absolute value of the correlation between the gene and the trait.
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Ch_22), sep="");
names(GSPvalue) = paste("p.GS.", names(Ch_22), sep="");

##############  Intramodular analysis: identifying genes with high GS and MM  #####################

####    Positive correlated   ######
module = "blue"                #0.65
module = "turquoise"           #0.79
module = "yellow"              #0.7  
module = "green"               #0.73

column = match(module, modNames);
moduleGenes = moduleColors == module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body Ch_22",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module);

#Intramodular analysis: identifying genes with high GS and MM
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
x <- abs(geneModuleMembership[moduleGenes, column])
y <- abs(geneTraitSignificance[moduleGenes, 1])
limit <- range(c(x,y)) 
r <- round(cor(x,y),2)
plot(    x, 
         y, 
         ylim =range(0.2,0.95), xlim =range(0.85, 0.95), 
         xlab =paste("Module Membership in", module, "module"), 
         ylab ="Gene significance for Ch22", col = module)
fit <-lm(y~x)
pintra = round(summary(fit)$coefficients[,4][[2]],2) #p-value
abline(fit, col='blue')
legend('topleft', col = "blue", lty = 1, box.lty = 1 ,legend = 'regression line', cex= 0.8) 
mtext(paste('correlation = ', r, ",", "P-value = ", pintra))
title(paste("Module membership vs. gene significance\n"))
text(y~x, labels = rownames(athalData3)[moduleColors=="blue"], data = athalData3, cex=0.8, font=4, pos = 2) #add label


# identifying genes with high GS and MM
intra_modular_analysis = data.frame(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(athalData3)[moduleColors=="yellow"] #only the x module
View(intra_modular_analysis)

# high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of the each driver gene in the x module to all other genes in the x module)
connectivity = intramodularConnectivity(adjacency, moduleColors)
connectivity = connectivity[colnames(athalData3)[moduleColors=="yellow"],] #only the x module
order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
connectivity = connectivity[order.kWithin,] #order rows following kWithin
connectivity = connectivity[1:5,] #top 5 genes that have a high connectivity to other genes in the x module
View(connectivity)
