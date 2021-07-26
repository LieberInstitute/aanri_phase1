#!/usr/bin/Rscript

library(tidyverse)
library(ggpubr)

get_pheno <- function(){
    fname = "../../../input/phenotypes/merged/_m/merged_phenotypes.csv"
    df = data.table::fread(fname) %>% select(-V1) %>%
        filter(Dx %in% c("Control"), Age > 17, Race %in% c("AA", "CAUC")) %>%
        mutate(Race = gsub("CAUC", "EA", Race))
    return(df)
}

memPHENO <- memoise::memoise(get_pheno)

get_cell_prop <- function(){
    load("../../_m/est_prop_Bisque.Rdata")
    cc = est_prop_bisque$caudate$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="Caudate")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="DLPFC")
    hh = est_prop_bisque$hippo$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="Hippocampus")
    gg = est_prop_bisque$dg$Est.prop.long %>%
        separate(sample, c("sample", "batch")) %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="Dentate Gyrus")
    df = bind_rows(cc, dd, hh, gg)
    return(df)
}

memPROP <- memoise::memoise(get_cell_prop)

get_cell_types <- function(){
    load("../../_m/est_prop_Bisque.Rdata")
    celltypes = unique(est_prop_bisque$caudate$Est.prop.long$cell_type)
    return(celltypes)
}

memCT <- memoise::memoise(get_cell_types)

save_img <- function(image, fn, w, h){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

pca_norm_data <- function(tissue){
    ## Load voom normalized data
    fname = paste0("../../../differential_analysis/", tissue,
                   "/_m/genes/voomSVA.RData")
    load(fname)
    ## Transpose expression
    norm_df = v$E %>% t
    ## Calculate PCA
    pca_df = prcomp(norm_df, center=TRUE)$x
    ## Convert to data frame
    norm_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        select(c(sample, PC1, PC2, PC3, PC4, PC5)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor)
    return(norm_dt)
}

pca_res_data <- function(tissue){
    ## Read in residualized data
    fname = paste0("../../../differential_analysis/", tissue,
                   "/_m/genes/residualized_expression.tsv")
    res_df = data.table::fread(fname) %>% column_to_rownames("V1") %>% t
    ## Calculate PCA
    pca_df = prcomp(res_df, center=TRUE)$x
    res_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        select(c(sample, PC1, PC2, PC3, PC4, PC5)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor)
    return(res_dt)
}

map_tissue <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus",
                "dentateGyrus"="Dentate Gyrus")[[tissue]])
}

for(tissue in c("caudate", "dlpfc", "hippocampus", "dentateGyrus")){
    dir.create(tissue)
    df <- memPROP() %>% filter(Tissue == map_tissue(tissue))
    for(ct in memCT()){
        ## Normalized
        sca = pca_norm_data(tissue) %>%
            inner_join(df, by="sample") %>% filter(cell_type == ct) %>%
            ggscatter(y="PC_values", x="Proportion", color="Race", palette="npg",
                      facet.by=c('PC'), ncol=5, add='reg.line', conf.int=TRUE,
                      cor.coef=TRUE, xlab=paste(ct, "Proportion"),
                      ylab="Normalized Expression",
                      panel.labs.font=list(face='bold', size = 14),
                      add.params=list(color="blue", fill="lightgray")) +
            font("xy.text", size=12) + font("xy.title", size=16, face="bold")
        save_img(sca, paste0(tissue,"/scatter_log2cpm_ancestry_5pcs_",ct),
                 w=18, h=6)
        ## Residualized
        sca = pca_res_data(tissue) %>%
            inner_join(df, by="sample") %>% filter(cell_type == ct) %>%
            ggscatter(y="PC_values", x="Proportion", color="Race", palette="npg",
                      facet.by=c('PC'), ncol=5, add='reg.line', conf.int=TRUE,
                      cor.coef=TRUE, xlab=paste(ct, "Proportion"),
                      ylab="Residualized Expression",
                      panel.labs.font=list(face='bold', size = 14),
                      add.params=list(color="blue", fill="lightgray")) +
            font("xy.text", size=12) + font("xy.title", size=16, face="bold")
        save_img(sca, paste0(tissue,"/scatter_resdf_ancestry_5pcs_",ct),
                 w=18, h=6)
    }
}

## Reproducibility Information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
