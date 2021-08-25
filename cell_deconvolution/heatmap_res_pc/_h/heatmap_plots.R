## Author: Louise Huuki
## Edited by KJ Benjamin

library(tidyverse)

## Functions
save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".pdf", ".png")){
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
        select(c(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(norm_dt)
}
memNORM <- memoise::memoise(pca_norm_data)

pca_res_data <- function(tissue){
    ## Read in residualized data
    fname = paste0("../../../differential_analysis/", tissue,
                   "/_m/genes/residualized_expression.tsv")
    res_df = data.table::fread(fname) %>% column_to_rownames("V1") %>% t
    ## Calculate PCA
    pca_df = prcomp(res_df, center=TRUE)$x
    res_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        select(c(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(res_dt)
}
memRES <- memoise::memoise(pca_res_data)

map_tissue <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus",
                "dentateGyrus"="Dentate Gyrus")[[tissue]])
}

get_pheno <- function(){
    fname = "../../../input/phenotypes/merged/_m/merged_phenotypes.csv"
    df = data.table::fread(fname) %>% select(-V1) %>%
        filter(Dx %in% c("Control"), Age > 17, Race %in% c("AA", "CAUC")) %>%
        mutate(Race = gsub("CAUC", "EA", Race))
    return(df)
}
memPHENO <- memoise::memoise(get_pheno)

get_cell_prop <- function(){
    ## Load Bisque Estimated Props
    load("../../_m/est_prop_Bisque.v2.Rdata")
    cc = est_prop_bisque$caudate$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="Caudate")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="DLPFC")
    hh = est_prop_bisque$hippo$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="Hippocampus")
    gg = est_prop_bisque$dg$Est.prop.long %>%
        separate(sample, c("sample", "batch")) %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="Dentate Gyrus")
    df = bind_rows(cc, dd, hh, gg)
    return(df)
}
memPROP <- memoise::memoise(get_cell_prop)

prep_data_prop <- function(tissue, func){
    df = memPROP() %>% filter(Tissue == map_tissue(tissue)) %>%
        mutate(RNum=sample)
    est_df <- inner_join(df, func(tissue), by="RNum") %>%
        select(RNum, cell_type, prop, PC, PC_values) %>%
        rename("Variable"="prop")
    est_df$PC <- factor(est_df$PC, levels = paste0('PC', 1:10))
    return(est_df)
}
memEST_PROP <- memoise::memoise(prep_data_prop)

prep_data <- function(covars, func, tissue){
    df = covars %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    est_df <- inner_join(df, func(tissue), by="RNum") %>%
        select(RNum, Covariate, Variable, PC, PC_values)
    est_df$PC <- factor(est_df$PC, levels = paste0('PC', 1:10))
    return(est_df)
}
memEST <- memoise::memoise(prep_data)

fit_model <- function(covars, fnc, tissue, celltype){
    ## Calculate p-values
    if(celltype){
        est_fit0 <- memEST_PROP(tissue, fnc) %>% group_by(cell_type, PC) %>%
            do(fitEST = broom::tidy(lm(Variable ~ PC_values, data=.)))
    } else {
        est_fit0 <- memEST(covars, fnc, tissue) %>% group_by(Covariate, PC) %>%
            do(fitEST = broom::tidy(lm(Variable ~ PC_values, data = .)))
    }
    est_fit = est_fit0 %>% unnest(fitEST) %>% filter(term == "PC_values") %>%
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf))
    print(est_fit %>% count(p.bonf.cat))
    return(est_fit)
}

tile_plot <- function(covars, fnc, tissue, fn, label, celltype=TRUE){
    ## Tile plot (heatmap)
    my_breaks <- c(0.05, 0.01, 0.005, 0)
    if(celltype){
        xlabel = "Cell Type"
        tile_plot <- fit_model(covars, fnc, tissue, celltype) %>%
            ggplot(aes(x = cell_type, y = PC, fill = log.p.bonf,
                       label = ifelse(p.bonf.sig,
                                      format(round(log.p.bonf,1), nsmall=1), "")))
        limits = c(0,55)
    } else {
        xlabel = "Covariate"
        tile_plot <- fit_model(covars, fnc, tissue, celltype) %>%
            ggplot(aes(x = Covariate, y = PC, fill = log.p.bonf,
                       label = ifelse(p.bonf.sig,
                                      format(round(log.p.bonf,1), nsmall=1), "")))
        limits = c(0,30)
    }
    tile_plot <- tile_plot + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        viridis::scale_color_viridis(option = "magma") +
        viridis::scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                                    direction=-1, limits=limits) +
        labs(x=xlabel, color="p-value Bonf\nsignificance",
             y=paste(label, "Expression (PCs)")) +
        ggpubr::theme_pubr(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    save_img(tile_plot, paste0(tissue,"/tilePlot_",fn))
}

#### Correlation with expression PCs ####
for(tissue in c("caudate", "dlpfc", "hippocampus", "dentateGyrus")){
    tname = tissue
    dir.create(tname)
    covarsCont = memPHENO() %>% select(-c("Region","BrNum","Sex","Race","Dx"))
    ## Normalized expression
    tile_plot(covarsCont, memNORM, tissue, "continous_norm", "Normalize", FALSE)
    tile_plot(covarsCont, memNORM, tissue, "celltypes_norm", "Normalize", TRUE)
    ## Residualized expression
    tile_plot(covarsCont, memRES, tissue, "continous_res","Residualized", FALSE)
    tile_plot(covarsCont, memRES, tissue, "celltypes_res","Residualized", TRUE)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
