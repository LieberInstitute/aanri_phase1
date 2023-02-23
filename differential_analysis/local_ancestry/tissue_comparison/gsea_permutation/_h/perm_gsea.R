## This script analyzsis gene set enrichment with permutation
suppressPackageStartupMessages({
    library(DOSE)
    library(here)
    library(dplyr)
    library(biomaRt)
    library(org.Hs.eg.db)
    library(clusterProfiler)
})

get_biomart <- function(){
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mapping <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),
                     mart = ensembl)
    return(mapping)
}
memBIOMART <- memoise::memoise(get_biomart)

get_mash_degs <- function(feature){
    fn <- paste0("../../summary_table/_m/",
                 "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    return(data.table::fread(fn) |>
           dplyr::filter(feature_type == feature) |>
           dplyr::mutate(length=abs(start - end)) |>
           dplyr::rename(beta=posterior_mean))
}

subset_tissue <- function(feature, tissue){
    return(get_mash_degs(feature) |>
           dplyr::filter(region == tissue) |>
           mutate(DE=ifelse(lfsr < 0.05, 1, 0)))
}

get_annotation <- function(feature){
    new_feature <- paste0(tolower(feature), "s")
    config      <- list("genes"="gene_annotation.tsv",
                        "transcripts"="tx_annotation.tsv",
                        "exons"="exon_annotation.tsv",
                        "junctions"="jxn_annotation.tsv")
    
    fn <- here("input/text_files_counts/_m/caudate", config[new_feature])
    dt <- data.table::fread(fn) |>
        mutate(ensembl_id=gsub("\\..*", "", gencodeID)) |>
        dplyr::rename("feature_id"="names") |>
        inner_join(memBIOMART(), by=c("ensembl_id"="ensembl_gene_id")) |>
        distinct(feature_id, .keep_all = TRUE) |>
        dplyr::select(feature_id, ensembl_id, entrezgene_id)
    return(dt)
}

add_entrez <- function(feature, tissue){
    return(subset_tissue(feature, tissue) |>
           inner_join(get_annotation(feature), by="feature_id"))
}

get_genes <- function(dt){
    return(dt |> dplyr::filter(DE == 1) |>
           dplyr::select(entrezgene_id, beta) |>
           mutate(abs_beta=abs(beta)) |> group_by(entrezgene_id) |>
           arrange(desc(abs_beta), .by_group=TRUE) |> slice(1) |>
           tidyr::drop_na() |> pull("entrezgene_id") |> as.character())
}

get_geneList <- function(dt){
    tmp <- dt |> na.exclude() |> dplyr::select(entrezgene_id, beta) |>
           mutate(abs_beta=abs(beta)) |> group_by(entrezgene_id) |>
           arrange(desc(abs_beta), .by_group=TRUE) |> slice(1)
    geneList <- tmp |> pull("beta")
    names(geneList) <- tmp$entrezgene_id
    return(geneList |> sort(decreasing=TRUE))
}

run_gseGO <- function(dt, ont, label){
                                        # GSEA with GO terms
    ego <- gseGO(geneList     = get_geneList(dt),
                 OrgDb        = org.Hs.eg.db,
                 ont          = ont,
                 keyType      = "ENTREZID",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 seed         = TRUE)
                                        # Save results
    if(dim(ego)[1] != 0){
        fn = paste0(label, "_", ont, ".tsv")
        ego <- setReadable(ego, "org.Hs.eg.db")
        ego |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_enichmentDGN <- function(dt, label){
    dgn <- enrichDGN(gene      = get_genes(dt),
                     universe  = names(get_geneList(dt)),
                     minGSSize = 5)
      # Save results
    if(dim(dgn)[1] != 0){
        fn = paste0(label, "_DGN.tsv")
        dgn <- setReadable(dgn, "org.Hs.eg.db")
        dgn |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_gseDGN <- function(dt, ont, label){
    dgn <- gseDGN(geneList     = get_geneList(dt),
                  minGSSize    = 5,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE,
                  seed         = TRUE)
      # Save results
    if(dim(dgn)[1] != 0){
        fn = paste0(label, "_", ont, ".tsv")
        dgn <- setReadable(dgn, "org.Hs.eg.db")
        dgn |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

##### MAIN
feature <- "Gene"
new_feature <- paste0(tolower(feature), "s")
dir.create(new_feature)
for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
    ##print(tissue)
    dt <- add_entrez(feature, tissue)
    label1 = paste0(new_feature,"/",gsub(" ", "_",tolower(tissue)),"_gsea")
    label2 = paste0(new_feature,"/",gsub(" ", "_",tolower(tissue)),"_enrich")
    for(ONT in c("BP", "MF", "CC")){
        run_gseGO(dt, ONT, label1)
    }
    run_gseDGN(dt, "DGN", label1)
    run_enichmentDGN(dt, label2)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
