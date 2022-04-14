## This script generates plots for region specific genes.
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggexport(p, filename=paste0(fn, ext), width=w, height=h)
    }
}

get_ancestry <- function(){
    file_name = paste0("../../../../input/ancestry_structure/structure.",
                       "out_ancestry_proportion_raceDemo_compare")
    return(data.table::fread(file_name))
}

discordant_genes <- function(){
    file_name = paste0("../../summary_table/region_specific/_m/",
                       "BrainSeq_ancestry_AA_shared_disconcordant.tsv")
    return(data.table::fread(file_name) %>% filter(!Concordant, Type == "Gene"))
}

get_biomart_df <- function(){
    return(data.table::fread("../_h/biomart.csv") %>%
           select(-c(V1, description)))
}
memMART <- memoise::memoise(get_biomart_df)

get_res_df <- function(){
    return(data.table::fread("residualized_expression.tsv") %>%
           filter(Geneid %in% discordant_genes()$Effect) %>%
           tibble::column_to_rownames("Geneid") %>% t %>% as.data.frame %>%
           tibble::rownames_to_column("RNum") %>%
           tidyr::pivot_longer(!RNum,names_to="Geneid",values_to="Expression"))
}
memRES <- memoise::memoise(get_res_df)

get_pheno <- function(){
    pheno_file <- "../../../../input/phenotypes/merged/_m/merged_phenotypes.csv"
    pheno_df   <- data.table::fread(pheno_file) %>%
        filter(Race == "AA", Age > 17, Dx == "Control") %>%
        inner_join(get_ancestry(), by=c("BrNum"="id"))
    return(pheno_df)
}
memPHENO <- memoise::memoise(get_pheno)

merge_data <- function(){
    return(memPHENO() %>% inner_join(memRES(), by="RNum") %>%
           mutate(ensembl_gene_id=gsub("\\..*", "", Geneid)) %>%
           inner_join(memMART(), by="ensembl_gene_id"))
}
memDF <- memoise::memoise(merge_data)

plot_scatter <- function(geneid){
    df  <- memDF() %>% filter(Geneid == geneid) %>%
        mutate(ID=ifelse(is.na(external_gene_name), Geneid, external_gene_name))
    sca <- ggscatter(df, x="Eur", y="Expression", color="Region", size=1,
                     alpha=0.6, add="reg.line", facet.by="external_gene_name",
                     palette="npg", panel.labs.font=list(face='bold'),
                     xlab="Genetic Ancestry (EA)",ylab="Residualized Expression",
                     ggtheme=theme_pubr(base_size=20, border=TRUE)) +
        font("xy.title", face="bold")
    fn = paste0("genetic_ancestry_discordant_", gsub("\\..*","",geneid))
    save_ggplots(fn, sca, 6, 6)
}

#### MAIN
gnames = memDF()$Geneid %>% unique
for(geneid in gnames){
    plot_scatter(geneid)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
