## This script generates meta plots form mash results
## for region specific.

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
})

combine_string <- function(glist){
    glist2 = c();
    for(xx in 1:length(strsplit(glist, "/"))){
        glist2 = c(glist2, strsplit(glist, "/")[[xx]])
    }
    return(unique(glist2))
}

get_immune_genes <- function(){
    fn <- here("differential_analysis/tissue_comparison",
               "gsea_permutation_pval1/_m",
               "immune_related_GSEA_results.tsv")
    dt <- data.table::fread(fn) %>% select(core_enrichment, Tissue)
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

gene_annotation <- function(){
                                        # Load significant DEGs
    fn <- here("differential_analysis/tissue_comparison/summary_table",
               "_m/BrainSeq_ancestry_4features_4regions.txt.gz")
    return(data.table::fread(fn) %>% filter(Type == "Gene") %>%
           select("Effect", "Symbol") %>% distinct)
}

get_eGenes_immune <- function(annot_df){
    fn <- "../../_m/BrainSeq_ancestry_4features_4regions.txt.gz"
    return(data.table::fread(fn) %>%
           filter(Feature == "Gene", gene_id %in% annot_df$Effect))
}

#### MAIN
immune_gnames <- get_immune_genes()
annot_df <- gene_annotation() %>% filter(Symbol %in% immune_gnames)

egene0 <- get_eGenes_immune(annot_df) %>% arrange(lfsr) %>%
    group_by(Tissue, gene_id) %>% slice(1) %>%
    select(Tissue, gene_id, variant_id, Symbol, lfsr, posterior_mean) %>%
    mutate(ID=ifelse(Symbol == "", gene_id, Symbol))

egene  <- egene0 %>% as.data.frame %>%
    select(Tissue, ID, posterior_mean) %>%
    tidyr::pivot_wider(names_from=Tissue, values_from=posterior_mean)

print(egene)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
