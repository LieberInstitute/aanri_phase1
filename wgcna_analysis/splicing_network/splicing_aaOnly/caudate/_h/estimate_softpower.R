## This script preforms WGCNA analysis
suppressMessages({
    library(dplyr)
    library(WGCNA)
})

PARAM_NETWORK_TYPE = 'signed'
options(stringsAsFactors = FALSE)
doParallel::registerDoParallel(cores=10)

filter_outliers <- function(expression, z_threshold = 2.5){
                                        # Input: an expression matrix
                                        # Output: an expression matrix with
                                        # outliers removed.
    sample_distance = dist(expression)
    dist_z = scale(colSums(as.matrix(sample_distance)))
    stopifnot(all(rownames(dist_z) == rownames(expression)))
    keepSamples = dist_z < z_threshold
    new_expression = expression[keepSamples,]
    new_expression
}

prepare_data <- function(){
                                        # Load sample data
    fn1 <- paste0("../../../../../differential_analysis/",
                  "caudate/_m/junctions/voomSVA.RData")
    load(fn1, verbose=TRUE)
    sample_table <- v$design %>% as.data.frame %>%
        select(-Intercept) %>% rename("Ancestry"="EA", "Sex"="Male")
                                        # Load residualized expression
    fn2 <- paste0("../../../../../differential_analysis/",
                  "caudate/_m/junctions/residualized_expression.tsv")
    vsd <- data.table::fread(fn2) %>% replace(is.na(.), "") %>%
        tibble::column_to_rownames("V1")
    print(dim(vsd))
                                        # Keep only the columns and rows that
                                        # are present in both the sample table
                                        # and vsd file
    samples      <- intersect(colnames(vsd), rownames(sample_table))
    vsd          <- vsd[,samples]
    sample_table <- sample_table[samples,]
                                        # WGCNA data import
    datExpr0 <- t(vsd)
                                        # Remove offending genes and samples
                                        # from the data
    gsg <- goodSamplesGenes(datExpr0, verbose = 3);
    if (!gsg$allOK){
        datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
    }
    datExpr <- datExpr0
                                        # Remove outliers
    datExpr <- filter_outliers(datExpr0, z_threshold = 2.5)
    rm(datExpr0)
                                        # Clean data
    samples      <- intersect(rownames(datExpr), rownames(sample_table))
    sample_table <- sample_table[samples,]
    datExpr      <- datExpr[samples,]
    print(dim(datExpr))
    save(datExpr, sample_table, file = '00.RData')
}


prepare_traits <- function(){
    load('00.RData', verbose=TRUE)
                                        # Associate traits with samples
    traitRows <- match(rownames(datExpr), rownames(sample_table))
    datTraits <- sample_table[traitRows,]
                                        # Diagnostic plot: Sample dendrogram and
                                        # trait heatmap
    pdf(file='sample_dendrogram_and_trait_heatmap.pdf',
        height=16, width = 22)
    sampleTree2 <- hclust(dist(datExpr), method = "average")
                                        # Convert traits to a color
                                        # representation: white means low, red
                                        # means high, grey means missing entry
    traitColors <- numbers2colors(traitRows, signed=FALSE);
                                        # Plot the sample dendrogram and the
                                        # colors underneath.
    plotDendroAndColors(sampleTree2, traitColors,
                        groupLabels="Avg. Counts",
                        main = "Sample dendrogram and trait heatmap",
                        cex.dendroLabels=0.7)
    dev.off()
    save(datExpr, sample_table, datTraits, file = "01.RData")
}

plot_power_parameter <- function(datExpr, plot_filename){
                                        # Choose a set of soft-thresholding
                                        # powers
    powers <- seq(1, 20, 1)
                                        # Call the network topology analysis
                                        # function
    sft = pickSoftThreshold(datExpr, networkType = PARAM_NETWORK_TYPE,
                            powerVector = powers, verbose = 5)
    print(sft$powerEstimate)
                                        # Plot the results:
    pdf(file=plot_filename)
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.5 , 2.2, 0.5), oma=c(0,0,2,0))
    cex1 = 0.7;
                                        # Scale-free topology fit index as a
                                        # function of the soft-thresholding power
    plot(sft$fitIndices[,1],
         -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1],
         -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="blue");
                                        # this line corresponds to using an R^2
                                        # cut-off of h
    abline(h=0.80,col="red")
                                        # Mean connectivity as a function of the
                                        # soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)", ylab="Mean Connectivity",
         type="n", main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5],
         labels=powers, cex=cex1, col="blue")
    #####
    plot(sft$fitIndices[,1], sft$fitIndices[,6],
         xlab="Soft Threshold (power)", ylab="Median Connectivity",
         type="n", main = paste("Median connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,6],
         labels=powers, cex=cex1, col="blue")
    #####
    plot(sft$fitIndices[,1], sft$fitIndices[,7],
         xlab="Soft Threshold (power)", ylab="Max Connectivity",
         type="n", main = paste("Max connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,7],
         labels=powers, cex=cex1, col="blue")
    dev.off()
}

figure_out_power_parameter <- function(){
    load(file = '01.RData', verbose=TRUE)
    plot_power_parameter(datExpr, 'power_parameter_selection.pdf')
}

#### MAIN
prepare_data()
                                        # Sample dendrogram and trait heatmap
prepare_traits()
                                        # Scale free topology model fit
figure_out_power_parameter()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
