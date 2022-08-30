## Expression correlation with cell types
library(tidyverse)

## Functions
save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".pdf", ".svg")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

map_tissue <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus",
                "dentateGyrus"="Dentate Gyrus")[[tissue]])
}

pca_norm_data <- function(region){
    baseloc <- paste0("/dcs04/lieber/statsgen/jbenjami/projects/",
                      "aanri_phase1/differential_analysis/")
    ## Load voom normalized data
    load(paste0(baseloc, region, "/_m/genes/voomSVA.RData"))
    ## Transpose expression
    norm_df = v$E %>% t
    ## Calculate PCA
    pca_df = prcomp(norm_df, center=TRUE)$x[, 1:20]
    ## Convert to data frame
    norm_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(norm_dt)
}
memNORM <- memoise::memoise(pca_norm_data)

pca_res_data <- function(region){
    baseloc <- paste0("/dcs04/lieber/statsgen/jbenjami/projects/",
                      "aanri_phase1/differential_analysis/")
    ## Read in residualized data
    fname = paste0(baseloc,region,"/_m/genes/residualized_expression.tsv")
    res_df = data.table::fread(fname) %>% column_to_rownames("V1") %>% t
    ## Calculate PCA
    pca_df = prcomp(res_df, center=TRUE)$x[, 1:20]
    res_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(res_dt)
}
memRES <- memoise::memoise(pca_res_data)

get_pheno <- function(region){
    baseloc <- paste0("/dcs04/lieber/statsgen/jbenjami/projects/",
                      "aanri_phase1/differential_analysis/")
    ## Load voom normalized data
    load(paste0(baseloc, region, "/_m/genes/voomSVA.RData"))
    dfx = v$design %>% as.data.frame %>% select(starts_with("qPC")) %>%
        rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    dfx$Covariate <- factor(dfx$Covariate, levels = paste0('qPC', 1:20))
    dfy = v$targets %>% as.data.frame %>%
        mutate(concordMapRate=mapply(function(r, n) {sum(r*n)/sum(n)},
                                     concordMapRate, numMapped)) %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor), as.numeric)) %>%
        mutate(across(where(is.logical), as.numeric)) %>%
        select(Afr,Sex,Age,RIN,mitoRate,rRNA_rate,totalAssignedGene,
               overallMapRate, concordMapRate)
    rownames(dfy) <- rownames(v$targets)
    dfy = dfy %>% rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    return(bind_rows(dfy, dfx) %>% mutate_if(is.character, as.factor))
}
memPHENO <- memoise::memoise(get_pheno)

get_cell_prop <- function(){
    ## Load Bisque Estimated Props
    load("../../_m/est_prop_Bisque.v2.Rdata")
    cc = est_prop_bisque$caudate$Est.prop.long %>%
        inner_join(memPHENO("caudate"), by=c("sample"="RNum")) %>%
        mutate(Tissue="Caudate")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>%
        inner_join(memPHENO("dlpfc"), by=c("sample"="RNum")) %>%
        mutate(Tissue="DLPFC")
    hh = est_prop_bisque$hippo$Est.prop.long %>%
        inner_join(memPHENO("hippocampus"), by=c("sample"="RNum")) %>%
        mutate(Tissue="Hippocampus")
    gg = est_prop_bisque$dg$Est.prop.long %>%
        separate(sample, c("sample", "batch")) %>%
        inner_join(memPHENO("dentateGyrus"), by=c("sample"="RNum")) %>%
        mutate(Tissue="Dentate Gyrus")
    df = bind_rows(cc, dd, hh, gg)
    return(df)
}
memPROP <- memoise::memoise(get_cell_prop)

prep_data_prop <- function(tissue, func){
    est_df <- memPROP() %>% mutate(RNum=sample) %>%
        inner_join(func(tissue), by="RNum") %>%
        select(RNum, cell_type, prop, PC, PC_values) %>%
        rename("Variable"="prop")
    est_df$PC <- factor(est_df$PC, levels = paste0('PC', 1:20))
    return(est_df)
}
memEST_PROP <- memoise::memoise(prep_data_prop)

fit_model <- function(){
    datalist1 <- list(); datalist2 <- list();
    fnc_lt <- list("Normalized"=memNORM, "Residualized"=memRES)
    for(fnc in c("Normalized", "Residualized")){
        for(tissue in c("caudate", "dentateGyrus", "dlpfc", "hippocampus")){
            est_fit <- memEST_PROP(tissue, fnc_lt[[fnc]]) %>%
                group_by(cell_type, PC) %>%
                do(fitEST = broom::tidy(lm(Variable~PC_values, data=.))) %>%
                unnest(fitEST) %>% filter(term != "(Intercept)") %>%
                mutate(p.bonf = p.adjust(p.value, "bonf"),
                       p.bonf.sig = p.bonf < 0.05,
                       p.bonf.cat = cut(p.bonf,
                                        breaks = c(1,0.05,0.01,0.005,0),
                                        labels = c("<= 0.005","<= 0.01",
                                                   "<= 0.05","> 0.05"),
                                        include.lowest = TRUE),
                       p.fdr = p.adjust(p.value, "fdr"),
                       log.p.bonf = -log10(p.bonf+10**(-300)))
            datalist1[[tissue]] <- est_fit %>%
                mutate(Tissue=map_tissue(tissue))
        }
        datalist2[[fnc]] <- bind_rows(datalist1) %>%
            mutate("Expression"=fnc)
    }
    return(bind_rows(datalist2))
}
memFIT <- memoise::memoise(fit_model)

tile_plot <- function(){
    my_breaks <- c(0.05, 0.01, 0.005, 0); limits = c(0,55)
    tile_plot <- memFIT() %>%
        ggplot(aes(x = cell_type, y = PC, fill = log.p.bonf,
                   label = ifelse(p.bonf.sig,
                                  format(round(log.p.bonf,1),
                                         nsmall=1),""))) +
        geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast=TRUE, aes(fontface="bold")) +
        facet_wrap("Tissue~Expression", ncol=4) +
        viridis::scale_color_viridis(option="magma") +
        viridis::scale_fill_viridis(name="-log10(p-value Bonf)",
                                    option="magma", direction=-1) +
        labs(x="Cell Type", y="Expression (PCs)") +
        ggpubr::theme_pubr(base_size=15) +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              panel.grid=element_blank(),
              axis.title=element_text(size=18, face="bold"),
              strip.text=element_text(face="bold", size=16))
    save_img(tile_plot, "tilePlot_celltypes_model_normalization", 16, 13)
}

#### MAIN ####
                                        # Correlation with expression PCs
memFIT() %>%
    data.table::fwrite("celltype_expression_normalization.tsv", sep='\t')

                                        # Tile plot (heatmap)
tile_plot()

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
