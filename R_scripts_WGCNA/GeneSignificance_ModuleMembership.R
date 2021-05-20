#################################################################################################
###
### Goal: Analysis of Gene Significance (GS) and Module Membership (MM)
###
### Made by: Cynthia Soto 
### Date: March 03, 2021
### Latest update: xxxx
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### DATA ASSUMPIONS:
### 1) Modules comes from   ** 17 SAMPLES OF A.THALIANA     I N F E C T E D     **  
### 2) WGCNA Std, perason corr, dissTOM Signed-Ntw

###################################################################################################
##
##                 Relating modules to external information and identifying importantgenes   
##                                   Gene Significance (GS) and ModuleMembership (MM)
##
###################################################################################################

#######################    BOTRYTIS CINEREA AT 24 hpi ############################################
## Define variable Bc24 containing the Bc24 column of datTrait

BotrytisC24 = as.data.frame(datTraits$Bc24)
sample = "BotrytisC24"
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

module = "blue"

## Jump to plot chunk 

#######################    BOTRYTIS CINEREA AT 12 hpi ############################################

BotrytisC12 = as.data.frame(datTraits$Bc12)
sample = "BotrytisC12"
names(BotrytisC12) = "BotrytisC12"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(athalData3, BotrytisC12, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(BotrytisC12), sep="")
names(GSPvalue) = paste("p.GS.", names(BotrytisC12), sep="")
## Intramodular analysis using the GS and MM values

module = "black"

## Jump to plot chunk 

#######################    BOTRYTIS CINEREA AT 18 hpi ############################################

BotrytisC18 = as.data.frame(datTraits$Bc18.1)
names(BotrytisC18) = "BotrytisC18"
sample = "BotrytisC18"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(athalData3, BotrytisC18, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(BotrytisC18), sep="")
names(GSPvalue) = paste("p.GS.", names(BotrytisC18), sep="")
## Intramodular analysis using the GS and MM values

module = "violet"

## Jump to plot chunk 

#######################    SCLEROTINIA 30 hpi ############################################

Sclerotinia30 = as.data.frame(datTraits$Ss30)
names(Sclerotinia30) = "Sclerotinia30"
sample = "Sclerotinia30"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(athalData3, Sclerotinia30, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Sclerotinia30), sep="")
names(GSPvalue) = paste("p.GS.", names(Sclerotinia30), sep="")
## Intramodular analysis using the GS and MM values

module = "blue"

## Jump to plot chunk 

#######################    COLLETOTRICHUM 40 hpi ############################################

Colletotrichum40 = as.data.frame(datTraits$Ch40)
names(Colletotrichum40) = "Colletotrichum40"
sample = "Colletotrichum40"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(athalData3, Colletotrichum40, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Colletotrichum40), sep="")
names(GSPvalue) = paste("p.GS.", names(Colletotrichum40), sep="")
## Intramodular analysis using the GS and MM values

module = "white"

##############################   Plot ###############################################################

column = match(module, modNames);
moduleGenes = moduleColors == module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));

title = paste("Module membership vs. gene significance for", module," module\n")
x_lab = paste("GS for A.thaliana with", sample)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = x_lab,
                   main = title,
                   cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = module)

