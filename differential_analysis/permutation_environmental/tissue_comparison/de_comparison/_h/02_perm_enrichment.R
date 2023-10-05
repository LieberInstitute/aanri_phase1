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

get_mash_degs <- function(){
    fn <- paste0("../../summary_table/_m/",
                 "BrainSeq_ancestry_binary_4features_4regions_allFeatures.txt.gz")
    return(data.table::fread(fn) |>
           dplyr::filter(feature_type == "Gene") |>
           dplyr::mutate(length=abs(start - end)) |>
           dplyr::rename(beta=posterior_mean))
}

subset_tissue <- function(tissue){
    return(get_mash_degs() |>
           dplyr::filter(region == tissue) |>
           mutate(DE=ifelse(lfsr < 0.05, 1, 0)))
}

get_annotation <- function(){
    new_feature <- "genes"
    config      <- list("genes"="gene_annotation.tsv",
                        "transcripts"="tx_annotation.tsv",
                        "exons"="exon_annotation.tsv",
                        "junctions"="jxn_annotation.tsv")
    
    fn <- here("input/text_files_counts/_m/caudate", config[new_feature])
    dt <- data.table::fread(fn) |>
        mutate(ensembl_id=gsub("\\..*", "", gencodeID)) |>
        dplyr::rename("feature_id"="names") |>
        inner_join(memBIOMART(), by=c("ensembl_id"="ensembl_gene_id"),
                   relationship = "many-to-many") |>
        distinct(feature_id, .keep_all = TRUE) |>
        dplyr::select(feature_id, ensembl_id, entrezgene_id)
    return(dt)
}

add_entrez <- function(tissue){
    return(subset_tissue(tissue) |>
           inner_join(get_annotation(), by="feature_id"))
}

get_comparison <- function(tissue){
    fn <- "../_m/BrainSeq_ancestry_DE_comparison.txt.gz"
    return(data.table::fread(fn) |>
           dplyr::filter(feature_type == "Gene", region == tissue))
}

get_genes <- function(dt, SHARED=FALSE){
    if(SHARED){
        genes <- dplyr::filter(get_comparison(tissue), ancestry == 1,
                               environ == 1)$feature_id
    } else {
        genes <- dplyr::filter(get_comparison(tissue), ancestry != 1,
                               environ == 1)$feature_id
    } 
    return(dt |> dplyr::filter(feature_id %in% genes) |>
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

run_enrichGO <- function(dt, ont, label, SHARED){
                                        # GSEA with GO terms
    ego <- enrichGO(gene     = get_genes(dt, SHARED),
                    universe = names(get_geneList(dt)),
                    OrgDb    = org.Hs.eg.db,
                    ont      = ont,
                    keyType  = "ENTREZID")
                                        # Save results
    if(dim(ego)[1] != 0){
        fn = ifelse(SHARED, paste0(label, "_", ont, ".shared.tsv"),
                    paste0(label, "_", ont, ".unique.tsv"))
        ego <- setReadable(ego, "org.Hs.eg.db")
        ego |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_enichmentDGN <- function(dt, label, SHARED){
    dgn <- enrichDGN(gene      = get_genes(dt, SHARED),
                     universe  = names(get_geneList(dt)),
                     minGSSize = 5)
                                        # Save results
    if(dim(dgn)[1] != 0){
        fn = ifelse(SHARED, paste0(label, "_DGN.shared.tsv"),
                   paste0(label, "_DGN.unique.tsv"))
        dgn <- setReadable(dgn, "org.Hs.eg.db")
        dgn |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

##### MAIN
outdir  <- "environmental"; dir.create(outdir)
for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
    print(tissue)
    dt <- add_entrez(tissue)
    label2 = paste0(outdir,"/",gsub(" ", "_",tolower(tissue)),"_enrich")
    for(shared in c(TRUE, FALSE)){
        for(ONT in c("BP", "MF", "CC")){
            run_enrichGO(dt, ONT, label2, shared)
        }
        run_enichmentDGN(dt, label2, shared)
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
