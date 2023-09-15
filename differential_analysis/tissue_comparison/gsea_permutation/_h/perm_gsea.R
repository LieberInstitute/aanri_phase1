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
    fn <- paste0("../../_m/", feature,"/lfsr_feature_4tissues.txt.gz")
    return(data.table::fread(fn) |>
           dplyr::rename("Dentate Gyrus"="Dentate.Gyrus"))
}

get_annotation <- function(feature){
    config = list("genes"="gene_annotation.tsv",
                  "transcripts"="tx_annotation.tsv",
                  "exons"="exon_annotation.tsv",
                  "junctions"="jxn_annotation.tsv")
    
    fn <- here("input/text_files_counts/_m/caudate/", config[feature])
    dt <- data.table::fread(fn) |>
        mutate(ensemblID=gsub("\\..*", "", gencodeID)) |>
        dplyr::rename("Effect"="names") |>
        inner_join(memBIOMART(), by=c("ensemblID"="ensembl_gene_id"),
                   relationship="many-to-many") |>
        distinct(Effect, .keep_all = TRUE)
    if(feature == "genes"){
        return(dt)
    } else {            
        return(dt |> mutate(Length=abs(start - end)))
    }
}

subset_tissue <- function(feature, tissue){
    return(get_mash_degs(feature) |>
           dplyr::select(c("Effect", all_of(tissue))) |>
           dplyr::rename("lfsr"=tissue) |>
           mutate(DE=ifelse(lfsr < 0.05, 1, 0)))
}

annot_effect_size <- function(feature, tissue){
    fn <- paste0("../../_m/", feature, "/posterior_mean_feature_4tissues.txt.gz")
    df <- data.table::fread(fn) |>
        dplyr::rename("Dentate Gyrus"="Dentate.Gyrus") |>
        dplyr::select(c("Effect", all_of(tissue))) |> dplyr::rename("beta"=tissue) |>
        inner_join(subset_tissue(feature, tissue), by="Effect")
    return(get_annotation(feature) |> inner_join(df, by="Effect"))
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
for(feature in c("genes", "transcripts", "exons", "junctions")){
    dir.create(feature)
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        ##print(tissue)
        dt <- annot_effect_size(feature, tissue)
        label1 = paste0(feature,"/",gsub(" ", "_",tolower(tissue)),"_gsea")
        label2 = paste0(feature,"/",gsub(" ", "_",tolower(tissue)),"_enrich")
        for(ONT in c("BP", "MF", "CC")){
            run_gseGO(dt, ONT, label1)
        }
        run_gseDGN(dt, "DGN", label1)
        run_enichmentDGN(dt, label2)
    }
}

                                        # Enrichment of all DEGs
feature = "genes"
for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
    print(tissue)
    dt <- annot_effect_size(feature, tissue)
    label2 = paste0(feature,"/",gsub(" ", "_",tolower(tissue)),"_enrich")
    for(ONT in c("BP", "MF", "CC")){
        print(ONT)
        run_enrichGO(dt, ONT, label2)
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
