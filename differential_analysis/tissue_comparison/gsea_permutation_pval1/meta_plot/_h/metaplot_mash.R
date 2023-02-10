## This script generates meta plots form mash results
## for region specific.

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(mashr)
})

get_biomart_df <- function(){
    return(data.table::fread(here("input/biomart/biomart.csv")) %>%
           select(-c(V1, description)) %>% tidyr::drop_na())
}

gene_annotation <- function(){
    fn = "../../../summary_table/_m/BrainSeq_ancestry_4features_4regions.txt.gz"
    return(data.table::fread(fn) %>% filter(Type == "Gene") %>%
           select("Effect", "Symbol") %>% distinct)
}

combine_string <- function(glist){
    glist2 = c();
    for(xx in 1:length(strsplit(glist, "/"))){
        glist2 = c(glist2, strsplit(glist, "/")[[xx]])
    }
    return(unique(glist2))
}

get_immune_genes <- function(){
    dt <- data.table::fread("../../_m/immune_related_GSEA_results.tsv") %>%
        select(core_enrichment, Tissue)
    cc <- dt %>% filter(Tissue == "Caudate") %>% pull(core_enrichment) %>%
        combine_string
    gg <- dt %>% filter(Tissue == "Dentate Gyrus") %>% pull(core_enrichment) %>%
        combine_string
    dd <- dt %>% filter(Tissue == "DLPFC") %>% pull(core_enrichment) %>%
        combine_string
    hh <- dt %>% filter(Tissue == "Hippocampus") %>% pull(core_enrichment) %>%
        combine_string
    return(intersect(intersect(intersect(cc, gg), dd), hh))
}

plot_meta <- function(geneid, gene_name){
    load("../../../_m/genes/mashr_meta_results.RData")
    labs <- c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
    libd.colors <- ggpubr::get_palette("npg", 4)
    plot.title <- gene_name
    fn <- paste("metaplot", plot.title, sep="_")
                                        # Print means
    print(gene_name)
    print(m$result$PosteriorMean[geneid,])
                                        # Plot meta and save
    R.devices::devEval(c("png", "pdf"), name=fn, {
        mash_plot_meta(m,get_significant_results(m)[geneid],ylab="",labels=labs,
                       colors=rmeta::meta.colors(
                                         box=as.character(libd.colors),
                                         lines="black", zero="black"))
        title(plot.title)
    })
}

#### MAIN
immune_gnames <- get_immune_genes()
annot_df <- gene_annotation() %>% filter(Symbol %in% immune_gnames)

for(ii in 1:dim(annot_df)[1]){
    plot_meta(annot_df[ii]$Effect, annot_df[ii]$Symbol)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
