## This script plots correlation between potential confounders.
suppressPackageStartupMessages({
    library(tidyverse)
})

## Functions
save_img <- function(image, fn, w=14, h=7){
    for(ext in c(".svg", ".pdf")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
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
    df = v$targets %>% as.data.frame %>%
        select(-c(BrNum, group, "norm.factors", "lib.size")) %>%
        mutate(concordMapRate=mapply(function(r, n) {sum(r*n)/sum(n)},
                                     concordMapRate, numMapped)) %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor), as.numeric)) %>%
        mutate(across(where(is.logical), as.numeric)) %>%
        select(Afr,Sex,Age,mitoRate,rRNA_rate,totalAssignedGene,
               overallMapRate, concordMapRate)
    rownames(df) <- rownames(v$targets)
    return(df)
}
memPHENO <- memoise::memoise(get_pheno)

prep_data <- function(covars){
    df = covars %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    return(df)
}
memEST <- memoise::memoise(prep_data)

merge_data <- function(covars, region){
    baseloc <- paste0("/dcs04/lieber/statsgen/jbenjami/projects/",
                      "aanri_phase1/input/phenotypes/_m/")
    qsv_lt  <- list("caudate"=paste0(baseloc, "qSV_caudate.csv"),
                    "dentateGyrus"=paste0(baseloc, "qSV_dg.csv"),
                    "dlpfc"=paste0(baseloc, "qSV_dlpfc.csv"),
                    "hippocampus"=paste0(baseloc, "qSV_hippo.csv"))
    qSV <- data.table::fread(qsv_lt[[region]]) %>% rename("RNum"="V1")
    qsv_df <- qSV %>% pivot_longer(!RNum, names_to="qSV", values_to="qSV_value")
    qsv_df$qSV <- factor(qsv_df$qSV, levels = paste0('PC', 1:15))
    df <- inner_join(memEST(covars), qsv_df, by="RNum")
    return(df)
}
memDF <- memoise::memoise(merge_data)

merge_expr <- function(covars, fnc, region){
    df <- inner_join(memEST(covars), fnc(region), by="RNum")
    df$PC <- factor(df$PC, levels = paste0('PC', 1:50))
    return(df)
}
memEXPR <- memoise::memoise(merge_expr)

merge_covars <- function(covars){
    df <- inner_join(memEST(covarsCont), memEST(covarsCont), by="RNum")
    return(df)
}
memCOVARS <- memoise::memoise(merge_covars)

fit_model <- function(covars, region, qsv, norm, identity){
    if(qsv){
        est_fit0 <- memDF(covars, region) %>% group_by(Covariate, qSV) %>%
            do(fitEST = broom::tidy(lm(Variable ~ qSV_value, data = .)))
    } else if(identity){
        est_fit0 <- memCOVARS(covars) %>% group_by(Covariate.x, Covariate.y) %>%
            do(fitEST = broom::tidy(lm(Variable.x ~ Variable.y, data = .)))
    } else {
        if(norm){
            est_fit0 <- memEXPR(covars, memNORM, region) %>%
                group_by(Covariate, PC) %>%
                do(fitEST = broom::tidy(lm(Variable ~ PC_values, data = .)))
        } else {
            est_fit0 <- memEXPR(covars, memRES, region) %>%
                group_by(Covariate, PC) %>%
                do(fitEST = broom::tidy(lm(Variable ~ PC_values, data = .)))
        }
    }
    ## Calculate p-values
    est_fit <- est_fit0 %>% unnest(fitEST) %>% filter(term != "(Intercept)") %>%
        mutate(p.bonf = p.adjust(p.value, "bonf"), p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf, breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf+10**(-300)))
    print(est_fit %>% count(p.bonf.cat))
    return(est_fit)
}
memFIT <- memoise::memoise(fit_model)

tile_plot <- function(covars, region, qsv=TRUE, norm=TRUE, identity=TRUE){
    ## Tile plot (heatmap)
    my_breaks <- c(0.05, 0.01, 0.005, 0)
    xlabel = "Covariate"
    if(qsv){
        ylabel = "qSV"; out = "qsv"
        tile_plot0 <- memFIT(covars, region, qsv, norm, identity) %>%
            ggplot(aes(x = Covariate, y = qSV, fill = log.p.bonf,
                       label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                      nsmall=1), "")))
        h = 7; w = 7; limits = c(0, 50)
    } else if(identity){
        ylabel = "Covariate"; out = "covars"
        tile_plot0 <- memFIT(covars, region, qsv, norm, identity) %>%
            mutate_if(is.character, as.factor) %>% rowwise() %>%
            mutate(pair=sort(c(Covariate.x,Covariate.y)) %>% paste(collapse=",")) %>%
            group_by(pair) %>% distinct(pair, .keep_all=TRUE) %>%
            ggplot(aes(x = Covariate.x, y = Covariate.y, fill = log.p.bonf,
                       label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                      nsmall=1), "")))
        h = 9; w = 9; limits = c(0, 100)
    } else {
        tile_plot0 <- memFIT(covars, region, qsv, norm, identity) %>%
            ggplot(aes(x = Covariate, y = PC, fill = log.p.bonf,
                       label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                      nsmall=1), "")))
        h = 10; w = 9; limits = c(0, 100)
        if(norm){
            ylabel = "Normalized Expression"; out = "norm"
        } else {
            ylabel = "Residualized Expression"; out = "res"
        }
    }
    tile_plot <- tile_plot0 + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE, aes(fontface="bold")) +
        viridis::scale_color_viridis(option = "magma") +
        viridis::scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                                    direction=-1, limits=limits) +
        labs(x=xlabel, y=ylabel) + ggpubr::theme_pubr(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid = element_blank())
    fn = paste0(region, "_tilePlot_",out,"_covariates")
    save_img(tile_plot, fn, w, h)
}

#### Correlation with expression PCs ####
for(region in c("caudate", "dentateGyrus", "dlpfc", "hippocampus")){
    outdir <- region
    ## Get covariates
    covarsCont = memPHENO(region) %>% rownames_to_column("RNum")
    ## Plotting
    tile_plot(covarsCont, region)
    tile_plot(covarsCont, region, FALSE, FALSE, TRUE)
    tile_plot(covarsCont, region, FALSE, TRUE, FALSE)
    tile_plot(covarsCont, region, FALSE, FALSE, FALSE)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
