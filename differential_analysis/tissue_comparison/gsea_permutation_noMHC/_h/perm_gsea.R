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

get_mhc_region <- function(mhc_region){
    mhc_lt <- list("mhc"="../../mhc_enrichment/_m/mhc_genes.csv",
                   "xmhc"="../../mhc_enrichment/extended_mhc/_m/xmhc_genes.csv",
                   "hla"="../../mhc_enrichment/hla_excluded/_m/hla_genes.csv")
    return(data.table::fread(mhc_lt[[mhc_region]]) |>
           dplyr::select(gene_id, gene_name) |>
           dplyr::mutate(ensembl_gene_id=gsub("\\..*", "", gene_id)))
}

remove_mhc_region <- function(mhc_region){
    return(memBIOMART() |>
           filter(!ensembl_gene_id %in% get_mhc_region(mhc_region)$ensembl_gene_id))
}

get_mash_degs <- function(feature){
    fn <- paste0("../../summary_table/_m/",
                 "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    return(data.table::fread(fn) |>
           dplyr::rename(beta=posterior_mean, feature_type=Type) |>
           dplyr::filter(feature_type == feature) |>
           dplyr::mutate(length=abs(start - end)))
}

subset_tissue <- function(feature, tissue){
    return(get_mash_degs(feature) |>
           dplyr::filter(Tissue == tissue) |>
           mutate(DE=ifelse(lfsr < 0.05, 1, 0)))
}

get_annotation <- function(feature, mhc_region){
    new_feature <- paste0(tolower(feature), "s")
    config      <- list("genes"="gene_annotation.tsv",
                        "transcripts"="tx_annotation.tsv",
                        "exons"="exon_annotation.tsv",
                        "junctions"="jxn_annotation.tsv")
    
    fn <- here("input/text_files_counts/_m/caudate", config[new_feature])
    dt <- data.table::fread(fn) |>
        mutate(ensembl_id=gsub("\\..*", "", gencodeID)) |>
        dplyr::rename("Effect"="names") |>
        inner_join(remove_mhc_region(mhc_region),
                   by=c("ensembl_id"="ensembl_gene_id"),
                   relationship="many-to-many") |>
        distinct(Effect, .keep_all = TRUE) |>
        dplyr::select(Effect, ensembl_id, entrezgene_id)
    return(dt)
}

add_entrez <- function(feature, tissue, mhc_region){
    return(subset_tissue(feature, tissue) |>
           inner_join(get_annotation(feature, mhc_region),
                      by="Effect"))
}

get_genes <- function(dt){
    return(dt |> dplyr::filter(DE == 1) |> dplyr::select(entrezgene_id, beta) |>
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

run_enrichGO <- function(dt, ont, label){
                                        # GSEA with GO terms
    ego <- enrichGO(gene     = get_genes(dt),
                    universe = names(get_geneList(dt)),
                    OrgDb    = org.Hs.eg.db,
                    ont      = ont,
                    keyType  = "ENTREZID")
                                        # Save results
    if(dim(ego)[1] != 0){
        fn = paste0(label, "_", ont, ".tsv")
        ego <- setReadable(ego, "org.Hs.eg.db")
        ego |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
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
set.seed(20230807)
for(mhc_region in c("hla", "mhc", "xmhc")){
    dir.create(mhc_region)
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        ##print(tissue)
        dt <- add_entrez("Gene", tissue, mhc_region)
        label1 = paste0(mhc_region,"/",gsub(" ", "_",tolower(tissue)),"_gsea")
        label2 = paste0(mhc_region,"/",gsub(" ", "_",tolower(tissue)),"_enrich")
        for(ONT in c("BP", "MF", "CC")){
            run_gseGO(dt, ONT, label1)
            run_enrichGO(dt, ONT, label2)
        }
        run_gseDGN(dt, "DGN", label1)
        run_enichmentDGN(dt, label2)
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
