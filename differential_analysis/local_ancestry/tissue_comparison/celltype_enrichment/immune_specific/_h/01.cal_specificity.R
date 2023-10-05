## Calculate PBMCs specificity.
## Code was adapted from
## https://github.com/jbryois/scRNA_disease/blob/master/Code_Paper/Code_Zeisel/get_Zeisel_Lvl4_input.md

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(Seurat)
    library(SingleCellExperiment)
})

#### Main

fn    <- here("input/public/randolph/inputs",
              "1_calculate_pseudobulk",
              "mergedAllCells_withCellTypeIdents_CLEAN.rds")
                                        # Load data
sce <- readRDS(fn)
sce <- as.SingleCellExperiment(sce)
keep_cells <- c("CD4_T", "CD8_T", "monocytes", "NK",
                "DC", "NKT", "B", "neutrophils")
sce <- sce[,(sce$celltype %in% keep_cells)]
colData(sce)$celltype <- droplevels(colData(sce)$celltype)
                                        # Aggregate dat
agg.sce <- scuttle::aggregateAcrossCells(sce, ids=sce$celltype,
                                         statistics="mean")
expr <- as.data.frame(counts(agg.sce)) |>
    tibble::rownames_to_column("gene_name")
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
expr <- expr |> group_by(gene_name) |>
    mutate(specificity=expr/sum(expr)) |>
    ungroup()
n_genes <- length(unique(expr$gene_name))
n_genes_to_keep <- round(n_genes * 0.1)
                                        # Extract markers based on specificity
dat <- expr |> filter(expr > 1) |> group_by(celltype) |>
    top_n(n=n_genes_to_keep, wt=specificity) |>
    ungroup() |> mutate(marker=1) |>
    select(gene_name, celltype, marker) |>
    tidyr::pivot_wider(names_from=celltype, values_from=marker) |>
    mutate_all(~tidyr::replace_na(., 0))
                                        # Save results
#data.table::fwrite(expr, "expression.pbmcs.tsv", sep="\t")
data.table::fwrite(dat, "randolph.single_cell.tsv", sep="\t")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
