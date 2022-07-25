## Gene ontology for RNA-seq
suppressPackageStartupMessages({
    library(dplyr)
    library(goseq)
    library(biomaRt)
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
    return(data.table::fread(fn) %>%
           dplyr::rename("Dentate Gyrus"="Dentate.Gyrus"))
}

get_annotation <- function(feature){
    config = list("genes"="gene_annotation.tsv",
                  "transcripts"="tx_annotation.tsv",
                  "exons"="exon_annotation.tsv",
                  "junctions"="jxn_annotation.tsv")
    
    fn <- paste0("/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/",
                 "input/text_files_counts/_m/caudate/", config[feature])
    dt <- data.table::fread(fn) %>%
        mutate(ensemblID=gsub("\\..*", "", gencodeID)) %>%
        dplyr::rename("Effect"="names") %>%
        inner_join(memBIOMART(), by=c("ensemblID"="ensembl_gene_id")) %>%
        distinct(Effect, .keep_all = TRUE)
    if(feature == "genes"){
        return(dt)
    } else {            
        return(dt %>% mutate(Length=abs(start - end)))
    }
}

subset_tissue <- function(feature, tissue){
    return(get_mash_degs(feature) %>%
           dplyr::select(c("Effect", all_of(tissue))) %>%
           dplyr::rename("lfsr"=tissue) %>%
           mutate(DE=ifelse(lfsr < 0.05, 1, 0)))
}

annot_effect_size <- function(feature, tissue){
    fn <- paste0("../../_m/", feature, "/posterior_mean_feature_4tissues.txt.gz")
    df <- data.table::fread(fn) %>%
        dplyr::rename("Dentate Gyrus"="Dentate.Gyrus") %>%
        dplyr::select(c("Effect", all_of(tissue))) %>% dplyr::rename("beta"=tissue) %>%
        inner_join(subset_tissue(feature, tissue), by="Effect")
    return(get_annotation(feature) %>% inner_join(df, by="Effect"))
}

get_goseq_genes <- function(dt, direction){
    if(tolower(direction) == "up"){
        genes <- dt %>% filter(beta > 0)
    } else if(tolower(direction) == "down"){
        genes <- dt %>% filter(beta < 0)
    } else {
        genes <- dt
    }
    return(genes %>% dplyr::select("ensemblID", "DE") %>%
           distinct(ensemblID, .keep_all=TRUE) %>%
           tibble::column_to_rownames("ensemblID") %>% pull("DE"))
}

get_gene_lengths <- function(dt, direction){
    if(tolower(direction) == "up"){
        genes <- dt %>% filter(beta > 0)
    } else if(tolower(direction) == "down"){
        genes <- dt %>% filter(beta < 0)
    } else {
        genes <- dt
    }
    return(genes %>% dplyr::select("ensemblID", "Length") %>%
           distinct(ensemblID, .keep_all=TRUE) %>%
           tibble::column_to_rownames("ensemblID") %>% pull("Length"))
}

run_goseq <- function(dt, direction, fn1, fn2){
    genes <- get_goseq_genes(dt, direction)
    gene_lengths <- get_gene_lengths(dt, direction)
    pdf(gsub(" ", "_", fn1))
    pwf <- nullp(genes, bias.data=gene_lengths, plot.fit=TRUE)
    dev.off()
                                        # Map GO
    new_fn1 <- gsub(" ", "_", paste(fn2,"kegg.txt", sep="_"))
    new_fn2 <- gsub(" ", "_", paste(fn2,"go.txt", sep="_")

    goseq(pwf, gene2cat=as.list(org.Hs.eg.db::org.Hs.egPATH)) %>%
        mutate(FDR_over=p.adjust(over_represented_pvalue, method="fdr"),
               FDR_under=p.adjust(under_represented_pvalue, method="fdr")) %>%
        data.table::fwrite(new_fn1, sep="\t")
    goseq(pwf, gene2cat=as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)) %>%
        mutate(FDR_over=p.adjust(over_represented_pvalue, method="fdr"),
               FDR_under=p.adjust(under_represented_pvalue, method="fdr")) %>%
        data.table::fwrite(new_fn2, sep="\t")
}

##### MAIN
for(feature in c("genes", "transcripts")){
    dir.create(feature)
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        ##print(tissue)
        dt <- annot_effect_size(feature, tissue)
        fn1 = paste0(feature,"/",tolower(tissue),"_pwf","_plotting.pdf")
        for(direction in c("ALL", "UP", "DOWN")){
            fn2 = paste0(feature, "/GOseq_analysis_mash_", tolower(tissue),
                         "_", tolower(direction))
            run_goseq(dt, direction, fn1, fn2)
        }
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
