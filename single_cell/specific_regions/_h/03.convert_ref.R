## Convert seurat object to h5ad

library(SingleCellExperiment)

load_data <- function(dataset){
    base_loc <- here::here("input/singlecell_datasets/_m",
                           "su_hpc_penn", dataset)
    meta  <- read.table(paste0(base_loc, "/meta.tsv"), header=TRUE,
                        sep="\t", as.is=TRUE, row.names=1)
    mat   <- data.table::fread(paste0(base_loc, "/exprMatrix.tsv.gz"))
    genes <- gsub(".+[|]", "", mat[,1][[1]])
    mat   <- data.frame(mat[,-1], row.names=genes)
    so    <- Seurat::CreateSeuratObject(counts=mat,
                                        project=dataset,
                                        meta.data=meta)
    return(Seurat::as.SingleCellExperiment(so))
}

### Main

                                        # Load R variable
for(dataset in c("microglia", "astrocyte", "oligo")){
    sce <- load_data(dataset)
    outfile <- paste0(dataset,".single_cell.h5ad")
    zellkonverter::writeH5AD(sce, file=outfile)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
