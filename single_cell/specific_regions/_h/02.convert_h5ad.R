## Convert seurat object to h5ad

library(SingleCellExperiment)

### Main

                                        # Load R variable
load("sce_aanri.brain_regions.rda")

                                        # Subset data
for(celltype in c("Micro", "Astro", "Oligo")){
    ct_sce  <- sce[, sce$new_celltypes == celltype]
                                        # Write as H5AD
    outfile <- paste0(tolower(celltype), ".aanri_brain_regions.h5ad")
    zellkonverter::writeH5AD(ct_sce, file=outfile)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
