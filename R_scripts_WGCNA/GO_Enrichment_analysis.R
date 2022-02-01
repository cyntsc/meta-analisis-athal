#####################################################################################
###
### Goal: WGCNA Std / Enrichement analysis 
###
### Made by: Cynthia Soto 
### Date: March 7, 2021
### Latest update: xxxx
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
### MODULES ARE FROM:
###                              ** 17     I N F E C T E D     ** 
###                                  SAMPLES OF A.THALIANA       
###
###
#####################################################################################

annot = read.csv(here("data", "geneGO_metadata.csv"), header=TRUE, sep='\t')
annot[1:5,]
dim(annot)
names(annot)

probes = names(athalData3);
probes2annot = match(probes, annot$Genes);

# Get the corresponding Locus Link IDs
allLLIDs = annot$Genes[probes2annot];

# Choose interesting 
intModules = c("green", "darkturquoise", "salmon")

for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

#####################################################################################
##
## enrichment analysis directly within R
##
#####################################################################################

## Biconductor packages GO.db, AnnotationDBI, and the appropriate organism-specific annotationpackage(s) need to be installed before running this code.
## The organism-specific packages have names of the formorg.Xx.eg.db, where Xx stands for organism code, for example, Mm for mouse, Hs for human, etc.  The only exceptionis yeast
## So this code needs the package org.Mm.eg.db

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "athaliana", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment

names(tab)

write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab


