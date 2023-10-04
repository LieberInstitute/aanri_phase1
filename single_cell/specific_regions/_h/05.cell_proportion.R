## Integrate datasets

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(ggpubr)
    library(Seurat)
    library(SingleCellExperiment)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

get_sc_data <- function(region){
    br_lt <- list("AMY"="SCE_AMY-n5_tran-etal.rda",
                  "DLPFC"="SCE_DLPFC-n3_tran-etal.rda",
                  "HPC"="SCE_HPC-n3_tran-etal.rda",
                  "NAc"="SCE_NAc-n8_tran-etal.rda",
                  "sACC"="SCE_sACC-n5_tran-etal.rda")
    fn <- here("input/singlecell_datasets/_m/tran",
               br_lt[[region]])
    load(fn)
    return(get(ls()[startsWith(ls(), "sce")]))
}

clean_data <- function(region){
                                        # Load data
    sce <- get_sc_data(region)
                                        # Remove "drop" cell types
    all_cells  <- as.character(unique(sce$cellType))
    keep_cells <- all_cells[!startsWith(all_cells, "drop")]
    sce <- sce[,(sce$cellType %in% keep_cells)]
    colData(sce)$cellType <- droplevels(colData(sce)$cellType)
                                        # Collapse similar cell types
    colData(sce)$new_celltypes <- stringr::str_split(colData(sce)$cellType,
                                                     "_", simplify=TRUE)[,1]
                                        # Remove prop data
    rowData(sce) <- rowData(sce)[!startsWith(names(rowData(sce)), "prop")]
                                        # Clean up
    reducedDims(sce) <- NULL; sizeFactors(sce) <- NULL
    metadata(sce)    <- list(NULL)
    sce$collapsedCluster <- NULL
                                        # Update object
    colData(sce)$region    <- region
    return(sce)
}

get_hvgs <- function(){
                                        # Load data
    sce.nac   <- clean_data("NAc")
    sce.hpc   <- clean_data("HPC")
    sce.dlpfc <- clean_data("DLPFC")
                                        # Biological components
    dec.nac   <- scran::modelGeneVar(sce.nac)
    dec.hpc   <- scran::modelGeneVar(sce.hpc)
    dec.dlpfc <- scran::modelGeneVar(sce.dlpfc)
    comb.dec  <- scran::combineVar(dec.nac, dec.hpc, dec.dlpfc)
    return(scran::getTopHVGs(comb.dec, n=5000))
}

load_data <- function(celltype){
                                        # Load data
    fn1   <- paste0(tolower(celltype), ".transfer_labels.h5ad")
    adata <- zellkonverter::readH5AD(fn1)
    fn2   <- paste0(tolower(celltype), ".aanri_brain_regions.h5ad")
    sce   <- zellkonverter::readH5AD(fn2)
                                        # Transfer labels
    colData(sce)$cell_type <- colData(adata)$C_scANVI
    colLabels(sce)         <- colData(sce)$cell_type
    return(sce)
}

calc_prop <- function(sce, transform="logit"){
    return(speckle::propeller(clusters=sce$cell_type,
                              sample=sce$sample_id,
                              group=sce$region,
                              transform=transform))
}

plot_prop <- function(sce, transform="logit"){
    props <- speckle::getTransformedProps(sce$cell_type,
                                          sce$sample_id,
                                          transform=transform)
    df  <- as.data.frame(props$TransformedProps) |>
        mutate(region = toupper(gsub("donor[0-9]_", "", sample)))
    bar <- ggbarplot(df, x="clusters", y="Freq", fill="region",
                     add="mean_se", palette="npg",
                     position=position_dodge(),
                     xlab="Celltypes", ylab="sin^-1(Proportion)",
                     ggtheme=theme_pubr(base_size=15))
    bxp <- ggboxplot(df, x="region", y="Freq", fill="region",
                     palette="npg", add="jitter", facet.by="clusters",
                     scales="free", panel.labs.font=list(face="bold"),
                     xlab="Celltypes", ylab="sin^-1(Proportion)",
                     ggtheme=theme_pubr(base_size=15, border=TRUE))
    return(list("BAR"=bar, "BOX"=bxp, "DATA"=df))
}

plot_clusters <- function(sce){
    set.seed(13)
    mnn.hold <- batchelor::fastMNN(sce, batch=sce$sample_id,
                                   d=50, correct.all=TRUE,
                                   subset.row=get_hvgs(),
                                   get.variance=TRUE,
                                   BSPARAM=BiocSingular::IrlbaParam())
    reducedDim(sce, "PCA_corrected_50") <- reducedDim(mnn.hold, "corrected")
    metadata(sce) <- metadata(mnn.hold)
    names(metadata(sce)) <- paste0(names(metadata(sce)),"_50")
    sce  <- scater::runTSNE(sce, dimred="PCA_corrected_50")
    sce  <- scater::runUMAP(sce, dimred="PCA_corrected_50")
    plt1 <- scater::plotReducedDim(sce, dimred="TSNE", colour_by="cell_type")
    plt2 <- scater::plotReducedDim(sce, dimred="UMAP", colour_by="cell_type")
    plt  <- plt1 + plt2
    return(plt)
}


#### Main

for(celltype in c("Micro", "Astro", "Oligo")){
    outfile0 <- paste0("propeller_proportion.",
                       tolower(celltype), ".anova.tsv")
    outfile1 <- paste0("transformed_proportions.",
                       tolower(celltype), ".tsv")
    sce      <- load_data(celltype)
    prop_df  <- calc_prop(sce, transform="asin")
    print(celltype)
    print(prop_df |> filter(FDR < 0.05))
    data.table::fwrite(prop_df, file=outfile0,
                       sep="\t", row.names=TRUE)
                                        # Generate plots
    plt      <- plot_clusters(sce)
    plot_lt  <- plot_prop(sce, "asin")
    data.table::fwrite(plot_lt[["DATA"]], file=outfile1, sep="\t")
                                        # Save plots
    outfile2 <- paste0(tolower(celltype), ".dim_reduction")
    outfile3 <- paste0(tolower(celltype), ".barplot")
    outfile4 <- paste0(tolower(celltype), ".boxplot")
    save_plot(plt, outfile2, 12, 6)
    save_plot(plot_lt[["BAR"]], outfile3, 10, 6)
    save_plot(plot_lt[["BOX"]], outfile4, 10, 6)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
