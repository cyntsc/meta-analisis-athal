#################################################################################################
###
### Goal: Enrichment analysis
###
### Made by: Cynthia Soto 
### Date: Feb 18, 2022
### Latest update: Feb xxx, 2022
###
### This is a PhD project associated at CICY Mx with the technical support (HPC) of 
###                                     Dr.Childs Lab at MSU.
###
#################################################################################################

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(athalData3, power = 9);
# Read in the annotation file
# annot = read.csv(file = "GeneAnnotation.csv");
annot = read.csv(here("data", "geneGO_metadata.csv"), header=TRUE, sep='\t')
annot[1:5,]
dim(annot)
names(annot)


# Select module
module = "darkmagenta";
# Select module probes
probes = names(athalData3)
inModule = (moduleColors==module);
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$Genes, annot$GO_term) )

nTop = 30;
IMConn = softConnectivity(athalData3[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$Genes, annot$GO_term) )
vis[1:10,1:5]

###############   Exporting Cytoscape

# Recalculate topological overlap if needed
# Select modules
modules = c("brown", "red");
modules = c("darkmagenta");

# Select module probes
probes = names(athalData3)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$GO_term[match(modProbes, annot$Genes)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
