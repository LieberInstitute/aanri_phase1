## This script compares popDE from limma::eBayes to mash model.
library(ggpubr)
library(magrittr)

get_eBayes <- function(){
    fn = "../../../deg_summary/_m/diffExpr_ancestry_full_4regions.tsv"
    return(data.table::fread(fn))
}
memBayes <- memoise::memoise(get_eBayes)

get_degs <- function(tissue, fdr=0.05){
    return(memBayes() %>%
           dplyr::filter(Tissue == tissue, adj.P.Val < fdr, Type == "Gene"))
}
memDEGs <- memoise::memoise(get_degs)

get_mash_model <- function(){
    load("../../_m/genes/mashr_meta_results.RData")
    return(list(posterior_mean = m$result$PosteriorMean %>% as.data.frame,
                posterior_sd = m$result$PosteriorSD %>% as.data.frame,
                lfsr = m$result$lfsr %>% as.data.frame))
}
memMODEL <- memoise::memoise(get_mash_model)

get_mash_tissue <- function(tissue, fdr=0.05){
    return(data.frame("Feature"=rownames(memMODEL()$lfsr),
                      "Mean"=memMODEL()$posterior_mean[, tissue],
                      "SD"=memMODEL()$posterior_sd[, tissue],
                      "lfsr"=memMODEL()$lfsr[, tissue]) %>%
           dplyr::filter(lfsr < fdr))
}
memMASH <- memoise::memoise(get_mash_tissue)

merge_data <- function(tissue){
    popDE_mash   <- memMASH(tissue) %>% dplyr::select(Feature, Mean)
    popDE_eBayes <- memDEGs(tissue, 1) %>% dplyr::select(Feature, logFC) %>%
        dplyr::filter(Feature %in% memMASH(tissue)$Feature)
    return(dplyr::inner_join(popDE_mash, popDE_eBayes, by="Feature"))
}
memDF <- memoise::memoise(merge_data)

save_ggplot <- function(fn, p, w=7, h=7){
    for(ext in c(".pdf", ".png")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

plot_effectsize <- function(tissue){
    sca <- ggscatter(memDF(tissue), x="Mean", y="logFC", add="reg.line",
                     add.params=list(color="blue", fill="lightgray"),
                     xlab="Effect size (mash)", ylab="Effect size (eBayes)",
                     conf.int=TRUE, cor.coef=TRUE,
                     cor.coeff.args=list(method="pearson", label.sep="\n"),
                     ggtheme=theme_pubr(base_size=15))
    return(sca)
}

##### MAIN
main <- function(){
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        fn  <- paste0("scatterplot_mash_eBayes_", gsub(" ","_",tolower(tissue)))
        sca <- plot_effectsize(tissue)
        save_ggplot(fn, sca)
    }
}

main()

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
