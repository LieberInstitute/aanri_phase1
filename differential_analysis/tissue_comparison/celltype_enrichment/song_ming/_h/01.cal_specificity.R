## Calculate specificity. Code was adapted from
## https://github.com/jbryois/scRNA_disease/blob/master/Code_Paper/Code_Zeisel/get_Zeisel_Lvl4_input.md

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(Seurat)
    library(SingleCellExperiment)
})

load_data <- function(dataset){
    base_loc <- here("input/singlecell_datasets/_m",
                     "su_hpc_penn", dataset)
    meta  <- read.table(paste0(base_loc, "/meta.tsv"), header=TRUE,
                       sep="\t", as.is=TRUE, row.names=1)
    mat   <- data.table::fread(paste0(base_loc, "/exprMatrix.tsv.gz"))
    genes <- gsub(".+[|]", "", mat[,1][[1]])
    mat   <- data.frame(mat[,-1], row.names=genes)
    so    <- CreateSeuratObject(counts=mat, project=dataset,
                                meta.data=meta)
    return(as.SingleCellExperiment(so))
}

agg_data <- function(dataset){
                                             # Load data
    sce <- load_data(dataset)
                                        # Aggregate dat
    agg.sce <- scuttle::aggregateAcrossCells(sce, ids=sce$subtype,
                                             statistics="mean")
    return(as.data.frame(counts(agg.sce)) |>
           tibble::rownames_to_column("gene_name"))
}

extract_spec <- function(dataset){
                                        # Aggregate data
    expr <- agg_data(dataset)
                                        # Filter not expressed genes
    not_expressed <- tidyr::pivot_longer(expr, !gene_name,
                                         names_to="celltype",
                                         values_to="expr") |>
        group_by(gene_name) |> summarize(total=sum(expr)) |>
        filter(total==0) |> pull(gene_name)
    expr <- filter(expr, !gene_name %in% not_expressed)
                                        # Scale data
    expr <- tidyr::pivot_longer(expr, !gene_name,
                                names_to="celltype",
                                values_to="expr") |>
        group_by(celltype) |>
        mutate(expr=expr*1e6/sum(expr))
                                        # Specificity calculation
    return(expr |> group_by(gene_name) |>
           mutate(specificity=expr/sum(expr)) |>
           ungroup())
}

extract_markers <- function(dataset){
    expr    <- extract_spec(dataset)
    n_genes <- length(unique(expr$gene_name))
    n_genes_to_keep <- round(n_genes * 0.1)
                                        # Extract markers based on specificity
    return(expr |> filter(expr > 1) |> group_by(celltype) |>
           top_n(n=n_genes_to_keep, wt=specificity) |>
           ungroup() |> mutate(marker=1) |>
           select(gene_name, celltype, marker) |>
           tidyr::pivot_wider(names_from=celltype, values_from=marker) |>
           mutate_all(~tidyr::replace_na(., 0)))
}
#### Main

for(dataset in c("microglia", "astrocyte", "oligo")){
    outfile <- paste0(dataset,".single_cell.tsv")
    extract_markers(dataset) |> data.table::fwrite(outfile, sep="\t")
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
