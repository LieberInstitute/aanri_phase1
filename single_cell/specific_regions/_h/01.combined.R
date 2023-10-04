## Integrate datasets
## modified from: https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_across-regions-analyses_step02_MNT.R

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(Seurat)
    library(SingleCellExperiment)
})

load_data <- function(region){
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
    sce <- load_data(region)
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

norm_data <- function(){
                                        # Load data
    sce.nac   <- clean_data("NAc")
    sce.hpc   <- clean_data("HPC")
    sce.dlpfc <- clean_data("DLPFC")
                                        # Combine datasets
    sce.aanriReg <- cbind(sce.nac, sce.hpc, sce.dlpfc)
    colData(sce.aanriReg)$sample_id <- paste0(sce.aanriReg$donor, "_",
                                                tolower(sce.aanriReg$region))
    sce.aanriReg$new_celltypes <- factor(sce.aanriReg$new_celltypes)
                                        # Re-calculate logcounts
    assay(sce.aanriReg, "logcounts") <- NULL
    sce.aanriReg <- batchelor::multiBatchNorm(sce.aanriReg,
                                                batch=sce.aanriReg$sample_id)
    return(sce.aanriReg)
}

#### Main
sce     <- norm_data()
prop_df <- speckle::propeller(clusters=sce$new_celltypes,
                              sample=sce$sample_id,
                              group=sce$region)
print(prop_df |> filter(FDR < 0.05))
data.table::fwrite(prop_df,
                   file="propeller_proportion.annriRegions.anova.tsv",
                   sep="\t", row.names=TRUE)
save(sce, file="sce_aanri.brain_regions.rda")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
