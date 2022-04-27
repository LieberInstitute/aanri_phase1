#########################################################
## This needs to be corrected across brain regions.    ##
## I will need to use dream instead here + cell props. ##
## Additional variables might also be needed.          ##
#########################################################

## Construct consensus models (pairwise) between brain regions
suppressMessages({
    library(dplyr)
    library(WGCNA)
})

PARAM_NETWORK_TYPE = 'signed'

get_cell_prop <- function(){
    celltype_file <- paste0("../../../../cell_deconvolution/_m/",
                            "est_prop_Bisque.v2.Rdata")
    load(celltype_file)
    cc = est_prop_bisque$caudate$Est.prop.long %>%
        mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>%
        mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    hh = est_prop_bisque$hippo$Est.prop.long %>%
        mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    gg = est_prop_bisque$dg$Est.prop.long %>%
        tidyr::separate(sample, c("sample", "batch")) %>%
        select(-batch) %>% mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    return(bind_rows(cc, dd, hh, gg) %>%
           tidyr::pivot_wider(names_from="cell_type", values_from="proportion"))
}
memPROP <- memoise::memoise(get_cell_prop)

get_ancestry <- function(){
    ancestry <- paste0("../../../../input/ancestry_structure/structure.",
                       "out_ancestry_proportion_raceDemo_compare")
    return(data.table::fread(ancestry) %>% select(-group))
}

get_phenotypes <- function(){
    fn <- "../../../../input/phenotypes/merged/_m/merged_phenotypes.csv"
    return(data.table::fread(fn) %>% select(-V1) %>%
           inner_join(get_ancestry(), by=c("BrNum"="id")) %>%
           inner_join(memPROP(), by="RNum"))
}
memPHENO <- memoise::memoise(get_phenotypes)

check_dup <- function(df){
    sample <- df %>% select_if(is.numeric)
    variables <- names(sample)
    return(cytominer::correlation_threshold(variables, sample, cutoff=0.95))
}

filter_outliers <- function(expression, z_threshold = 2.5){
     # Input: an expression matrix
     # Output: an expression matrix with outliers removed
     # Remove samples with z normalized total distance from other samples > z_threshold
     sample_distance = dist(expression)
     dist_z = scale(colSums(as.matrix(sample_distance)))
     stopifnot(all(rownames(dist_z) == rownames(expression)))
     keepSamples = dist_z < z_threshold
     new_expression = expression[keepSamples,]
     new_expression
}

prepare_data <- function(setLabels){
    ## Load sample data
    sample_table <- memPHENO() %>% select(-Afr) %>%
        filter(Age > 17, Race %in% c("AA", "CAUC")) %>%
        mutate_if(is.numeric, scales::rescale) %>%
        tibble::column_to_rownames("RNum")
    if(length(check_dup(sample_table)) != 0){
        sample_table <- sample_table %>% select(-check_dup(sample_table))
    }

    ## Load residualized expression
    fn  <- paste0("../../../../internal_replication/tissue_comparison/",
                  "residualized_expression/_m/genes/residualized_expression.tsv")
    vsd <- data.table::fread(fn) %>% replace(is.na(.), "") %>%
        tibble::column_to_rownames("Geneid")
    print(dim(vsd))

    ## Keep only the columns and rows that are present in
    ## both the sample table and vsd file
                                        # Caudate
    samples_cc <- intersect(colnames(vsd),
                            rownames(sample_table %>% filter(Region == "Caudate")))
    vsd_cc = vsd[,samples_cc]
                                        # Dentate Gyrus
    samples_gg <- intersect(colnames(vsd),
                            rownames(sample_table %>% filter(Region == "DentateGyrus")))
    vsd_gg = vsd[,samples_gg]
                                        # DLPFC
    samples_dd <- intersect(colnames(vsd),
                            rownames(sample_table %>% filter(Region == "DLPFC")))
    vsd_dd = vsd[,samples_dd]
                                        # Hippocampus
    samples_hh <- intersect(colnames(vsd),
                            rownames(sample_table %>% filter(Region == "HIPPO")))
    vsd_hh = vsd[,samples_hh]

    ## WGCNA data import
    nSets = 4
    shortLabels = c("Caudate", "Dentate_Gyrus", "DLPFC", "Hippocampus")
    multiExpr0 = vector(mode="list", length=nSets)
                                        # Caudate
    multiExpr0[[1]] = list(data=as.data.frame(t(vsd_cc)))
    names(multiExpr0[[1]]$data) = rownames(vsd_cc)
    rownames(multiExpr0[[1]]$data) = colnames(vsd_cc)
                                        # Dentate Gyrus
    multiExpr0[[2]] = list(data=as.data.frame(t(vsd_gg)))
    names(multiExpr0[[2]]$data) = rownames(vsd_gg)
    rownames(multiExpr0[[2]]$data) = colnames(vsd_gg)
                                        # DLPFC
    multiExpr0[[3]] = list(data=as.data.frame(t(vsd_dd)))
    names(multiExpr0[[3]]$data) = rownames(vsd_dd)
    rownames(multiExpr0[[3]]$data) = colnames(vsd_dd)
                                        # Hippocampus
    multiExpr0[[4]] = list(data=as.data.frame(t(vsd_hh)))
    names(multiExpr0[[4]]$data) = rownames(vsd_hh)
    rownames(multiExpr0[[4]]$data) = colnames(vsd_hh)
                                        # Check data
    exprSize = checkSets(multiExpr0)
    print(exprSize)
    # Remove offending genes and samples from the data
    gsg = goodSamplesGenesMS(multiExpr0, verbose = 3);
    if (!gsg$allOK){
        for(set in 1:exprSize$nSets){
            multiExpr0[[set]]$data = multiExpr0[[set]]$data[gsg$goodSamples, gsg$goodGenes]
        }
    }
    # Secondary sample filtering
    for(set in 1:exprSize$nSets){
        multiExpr0[[set]]$data = filter_outliers(multiExpr0[[set]]$data, 2.5)
    }
    multiExpr <- multiExpr0
    exprSize = checkSets(multiExpr)
    samples_cc = intersect(rownames(multiExpr[[1]]$data), samples_cc)
    samples_gg = intersect(rownames(multiExpr[[2]]$data), samples_gg)
    samples_dd = intersect(rownames(multiExpr[[3]]$data), samples_dd)
    samples_hh = intersect(rownames(multiExpr[[4]]$data), samples_hh)
    samples = c(samples_cc, samples_gg, samples_dd, samples_hh)
    sample_table = sample_table[samples,]
    save(multiExpr, exprSize, sample_table, shortLabels, file = '00.RData')
}

plot_sample_clustering <- function(setLabels){
    lnames = load('00.RData')
    sampleTrees = list()
    for(set in 1:exprSize$nSets){
        sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method="average")
    }
    pdf(file='sample_clustering.pdf', height=12, width=12)
    par(mfrow=c(2,1))
    par(mar=c(0,4,2,0))
    for(set in 1:exprSize$nSets){
        plot(sampleTrees[[set]],
             main=paste("Sample clustering on all genes in ", setLabels[set]),
             xlab="", sub="", cex=0.7)
    }
    dev.off()
}

prepare_traits <- function(){
    lnames = load('00.RData')
    Traits <- vector(mode="list", length=exprSize$nSets)
    # Associate traits with samples
    for(set in 1:exprSize$nSets){
        setSamples = rownames(multiExpr[[set]]$data)
        traitRows = match(setSamples, rownames(sample_table))
        Traits[[set]] = list(data=sample_table[traitRows, c(-1)])
        rownames(Traits[[set]]$data) = rownames(sample_table[traitRows, ])
    }
    nGenes = exprSize$nGenes
    nSamples = exprSize$nSamples
    save(multiExpr, exprSize, sample_table, shortLabels,
         Traits, nGenes, nSamples, file = "01.RData")
}

plot_power_parameter <- function(nSets, multiExpr, RsquaredCut = 0.85){
    # Choose a set of soft-thresholding powers
    powers = seq(from = 4, to=20, by=1)
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets)
    softPowerTables = vector(mode = "list", length = nSets)
    # Call the network topology analysis function for each set in turn
    for (set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data,
                                                           powerVector=powers, verbose = 2,
                                                           networkType=PARAM_NETWORK_TYPE)[[2]])
        # Calculated softpower from fitted values
        cond = powerTables[[set]]$data$`SFT.R.sq` > RsquaredCut
        softPowerTables[[set]] = min(powerTables[[set]]$data[cond,"Power"])
    }
    ##softpower = max(unlist(softPowerTables))
    softpower = 12 # based on no convergence for one brain region (HIPPO??)
    print(softpower)
    # Plot the results:
    colors = c("black", "red")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity",
                 "Median connectivity", "Max connectivity")
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4)
    for (set in 1:nSets){
        for (col in 1:length(plotCols)){
            ylim[1, col] = min(ylim[1, col],
                               powerTables[[set]]$data[, plotCols[col]],
                               na.rm = TRUE)
            ylim[2, col] = max(ylim[2, col],
                               powerTables[[set]]$data[, plotCols[col]],
                               na.rm = TRUE)
        }
    }
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    pdf(file = "power_parameter_selection.pdf", wi = 8, he = 6)
    sizeGrWindow(8, 6)
    par(mfcol = c(2,2))
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7
    for (col in 1:length(plotCols)) for (set in 1:nSets){
        if (set==1){
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
                 main = colNames[col])
            addGrid()
        }
        if (col==1){
            text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 labels=powers,cex=cex1,col=colors[set])
        } else {
            text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                 labels=powers,cex=cex1,col=colors[set])
        }
        if (col==1){
            legend("bottomright", legend = setLabels, col = colors, pch = 20)
        } else {
            legend("topright", legend = setLabels, col = colors, pch = 20)
        }
    }
    dev.off()
    return(softpower)
}

figure_out_power_parameter <- function(){
    lnames = load('01.RData')
    nSets = exprSize$nSets
    softpower <- plot_power_parameter(nSets, multiExpr, 0.85)
    return(softpower)
}

construct_network <- function(softPower){
    enableWGCNAThreads()
    lnames = load("01.RData")
    ## softPower value from previous plot power_parameter_selection.pdf
    cor <- WGCNA::cor
    net = blockwiseConsensusModules(multiExpr, maxBlockSize=30000,
                                    power=softPower, minModuleSize=30,
                                    deepSplit=2, pamRespectsDendro=FALSE,
                                    mergeCutHeight=0.25, numericLabels=TRUE,
                                    minKMEtoStay=0, corType="bicor",
                                    saveTOMFileBase="TOM", saveTOMs=TRUE,
                                    networkType=PARAM_NETWORK_TYPE,
                                    TOMType=PARAM_NETWORK_TYPE, verbose=3)
    consMEs = net$multiMEs
    moduleLabels = net$colors
    moduleColors = labels2colors(moduleLabels)
    consTree = net$dendrograms[[1]]
    save(net, consMEs, moduleLabels, moduleColors, consTree, file="02.RData")
}

plot_cluster_dendrogram <- function(){
    lnames = load("02.RData")
    pdf(file = "consensus_dendrogram.pdf", wi = 8, he = 6)
    sizeGrWindow(8,6)
    plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels=FALSE,
                        hang=0.03, addGuide=TRUE, guideHang=0.05,
                        main="Consensus gene dendrogram and module colors")
    dev.off()
}

consensus_eigengene_network <- function(){
    lnames = load(file = "01.RData")
    lnames = load(file = "02.RData")
    nSets = exprSize$nSets
    # Create a variable weight that will hold just the body weight of mice in both sets
    region = vector(mode = "list", length = nSets);
    for (set in 1:nSets){
        region[[set]] = list(data = as.data.frame(Traits[[set]]$data$Region))
        names(region[[set]]$data) = "Region"
    }
    # Recalculate consMEs to give them color names
    consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)
    # Plot eigengene network
    pdf(file = "eigengene_networks.pdf", width=8, height=10)
    sizeGrWindow(8,10)
    par(cex = 0.9)
    plotEigengeneNetworks(consMEsC, setLabels, marDendro=c(0,2,2,1),
                          marHeatmap=c(3,3,2,1), xLabelsAngle=0,
                          zlimPreservation=c(0.5, 1))
    dev.off()
    # We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
    MET = consensusOrderMEs(addTraitToMEs(consMEsC, region))
    # Plot eigengene network
    pdf(file = "eigengene_networks_region.pdf", width=8, height=10)
    sizeGrWindow(8,10)
    par(cex = 0.9)
    plotEigengeneNetworks(MET, setLabels, marDendro=c(0,2,2,1),
                          marHeatmap=c(3,3,2,1), xLabelsAngle=0,
                          zlimPreservation=c(0.5, 1))
    dev.off()
    save(MET, consMEsC, region, file="03.RData")
}

export_eigengene_tables <- function(){
    lnames = load(file = "01.RData")
    lnames = load(file = "02.RData")
    lnames = load(file = "03.RData")
    nSets = exprSize$nSets
    ## Export eigengene tables
    for(set in 1:nSets){
        write.csv(consMEsC[[set]]$data,
                  paste0('eigengenes_',shortLabels[[set]],'.csv'))
    }
    # Write modules
    modules = data.frame(row.names=colnames(multiExpr[[1]]$data),
                         module=moduleColors)
    write.csv(modules, 'modules.csv')
}

#### MAIN
setLabels = c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
prepare_data(setLabels)
plot_sample_clustering(setLabels)
prepare_traits()
softpower <- figure_out_power_parameter()
construct_network(softpower)
plot_cluster_dendrogram()
consensus_eigengene_network()
export_eigengene_tables()


#### Reproducibility Information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
