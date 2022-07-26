## This script preforms WGCNA analysis
suppressMessages({
    library(dplyr)
    library(WGCNA)
})

PARAM_NETWORK_TYPE = 'signed'
options(stringsAsFactors = FALSE)
doParallel::registerDoParallel(cores=15)

construct_network <- function(softPower){
    load(file = "01.RData", verbose=TRUE)
    enableWGCNAThreads(nThreads = 15)
    cor <- WGCNA::cor
    net = blockwiseModules(datExpr, mergeCutHeight = 0.2,
                           power = softPower, minModuleSize = 30,
                           networkType = PARAM_NETWORK_TYPE,
                           TOMType = PARAM_NETWORK_TYPE,
                           numericLabels = TRUE, corType = "bicor",
                           saveTOMs = TRUE, saveTOMFileBase = "TOM",
                           verbose = 3, maxBlockSize=28750)
    moduleLabels = net$colors
    moduleColors = labels2colors(net$colors)
    MEs = net$MEs;
    geneTree = net$dendrograms[[1]];
    save(net, MEs, moduleLabels, moduleColors, geneTree,
         softPower, file="02.RData")
}

plot_cluster_dendrogram <- function(){
    load(file = "02.RData", verbose=TRUE)
    pdf(file="cluster_dendrogram.pdf", height=16, width=22)
    mergedColors = labels2colors(net$colors)
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                        "Module Colors", dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05, cex.dendroLabels=0.3)
    dev.off()
}

correlate_with_traits <- function(){
    load(file = "01.RData", verbose=TRUE)
    load(file = "02.RData", verbose=TRUE)
                                        # Define numbers of genes and samples
    nGenes = ncol(datExpr); nSamples = nrow(datExpr);
                                        # Recalculate MEs with color labels
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, datTraits, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
                                        # Plot
    pdf(file="module_trait_relationships.pdf", height=22,width = 26)
                                        # Will display correlations and their
                                        # p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3));
                                        # Display the correlation values within
                                        # a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),
                   yLabels = names(MEs), ySymbols = names(MEs),
                   colorLabels = FALSE, naColor = "grey",
                   colors = blueWhiteRed(50), textMatrix = textMatrix,
                   setStdMargins = FALSE, cex.text = 0.9, zlim = c(-1,1),
                   main = paste("Module kME-Trait Correlation"))
    dev.off()
}

export_eigengene_tables <- function(){
    load(file = "01.RData", verbose=TRUE)
    load(file = "02.RData", verbose=TRUE)
                                        # Define numbers of genes and samples
    nGenes   <- ncol(datExpr)
    nSamples <- nrow(datExpr)
                                        # Recalculate MEs with color labels
    MEs0           <- moduleEigengenes(datExpr, moduleColors)$eigengenes
    rownames(MEs0) <- rownames(datExpr)
    write.csv(MEs0, 'eigengenes.csv')
                                        # Write modules
    modules <- data.frame(row.names=colnames(datExpr), module=moduleColors)
    write.csv(modules, 'modules.csv')
    save(datExpr,softPower,moduleColors, file = "cytoscapenetwork.Rdata")
}

#### MAIN
                                        # Softpower based on Dentate Gyrus
construct_network(softPower=15)
                                        # TOM dendrogram
plot_cluster_dendrogram()
                                        # Module eigenvalue correlation
correlate_with_traits()
export_eigengene_tables()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
