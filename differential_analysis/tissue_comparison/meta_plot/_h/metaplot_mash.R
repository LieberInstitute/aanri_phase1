## This script generates meta plots form mash results
## for region specific.

suppressPackageStartupMessages({
    library(dplyr)
    library(mashr)
})

get_biomart_df <- function(){
    return(data.table::fread("../../../../input/biomart/biomart.csv") %>%
           select(-c(V1, description, entrezgene)))
}
memMART <- memoise::memoise(get_biomart_df)

get_region_specific <- function(){
    file_name = paste0("../../summary_table/region_specific/_m/",
                       "BrainSeq_ancestry_AA_region_specific.tsv")
    return(data.table::fread(file_name) %>% filter(Type == "Gene"))
}

get_genes_by_effect <- function(){
    return(get_region_specific()%>% arrange(Tissue,desc(abs(posterior_mean)))%>%
           group_by(Tissue) %>% slice(1) %>%
           mutate(ensembl_gene_id=strsplit(gencodeID, "[.]")[[1]][1]) %>%
           inner_join(memMART(), by="ensembl_gene_id"))
}

get_genes_by_lfsr_N_effect <- function(){
    return(get_region_specific() %>% arrange(Tissue, lfsr) %>%
           group_by(Tissue) %>% slice(1:5) %>%
           arrange(Tissue, desc(abs(posterior_mean))) %>%
           group_by(Tissue) %>% slice(1) %>%
           mutate(ensembl_gene_id=strsplit(gencodeID, "[.]")[[1]][1]) %>%
           inner_join(memMART(), by="ensembl_gene_id"))
}

plot_meta <- function(geneid, gene_name){
    load("../../_m/genes/mashr_meta_results.RData")
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
geneids1 <- get_genes_by_effect()$Effect
gnames1  <- get_genes_by_effect()$external_gene_name
geneids2 <- get_genes_by_lfsr_N_effect()$Effect
gnames2  <- get_genes_by_lfsr_N_effect()$external_gene_name

for(ii in seq(1,4)){
    plot_meta(geneids1[[ii]], gnames1[[ii]])
    plot_meta(geneids2[[ii]], gnames2[[ii]])
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
