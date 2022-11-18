## mash modeling with common baseline at mean
##
## Author: Kynon J Benjamin
##
#############################################

suppressPackageStartupMessages({
    library(mashr)
    library(dplyr)
    library(ggpubr)
    library(argparse)
})

save_img <- function(image, fn, w, h){
    for(ext in c(".svg", ".pdf")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_bhat <- function(feature){
    fn = paste0(feature, "/bhat_de_4tissues.tsv")
    return(data.table::fread(fn) %>% mutate(effect=gene_id) %>%
           distinct(effect, .keep_all=TRUE))
}

get_shat <- function(feature){
    fn = paste0(feature, "/shat_de_4tissues.tsv")
    return(data.table::fread(fn) %>% mutate(effect=gene_id) %>%
           distinct(effect, .keep_all=TRUE))
}

plot_mixture_prop <- function(m, feature){
    fn = paste0(feature,"/barplot_estimated_pi")
    df = get_estimated_pi(m) %>% as.data.frame %>%
        tibble::rownames_to_column("Model")
    colnames(df)[2] <- "Estimated pi"
    brp = ggbarplot(df, x="Model", y="Estimated pi", fill="gray",
                    ggtheme=theme_pubr(base_size=20), xlab="",
                    label=TRUE, label.pos="out", lab.nb.digits=2) +
        font("y.title", face="bold") + rotate_x_text(45)
    save_img(brp, fn, 7, 7)
}

run_mashr <- function(feature){
    ## Load prepared data
    bhat <- get_bhat(feature) %>%
        tibble::column_to_rownames("effect") %>%
        select(-gene_id) %>% as.matrix
    shat <- get_shat(feature) %>%
        tibble::column_to_rownames("effect") %>%
        select(-gene_id) %>% as.matrix
    ## Prepare for mashr
    data.init = mash_set_data(bhat, shat)
    Vhat = estimate_null_correlation_simple(data.init)
    ## Update mash
    data = mash_update_data(data.init, V=Vhat)
    rm(data.init)
    ## Calculate data driven covariances
    U.c = cov_canonical(data)
                                        # Condition by condition for strong set
    m.1by1 = mash_1by1(data)
    strong = get_significant_results(m.1by1)
    rm(m.1by1)
    U.pca = cov_pca(data, 4, subset=strong)
    U.ed = cov_ed(data, U.pca, subset=strong)
    rm(strong, U.pca)
    ## Fit mash model and compute posterior summaries
    m = mash(data, Ulist = c(U.ed,U.c), algorithm.version="Rcpp")
    rm(U.ed, U.c)
    ## Examine pairwise sharing
    print(get_pairwise_sharing(m))
    ## Save significant results
    print(get_significant_results(m) %>% length)
    save(m, data, file=paste0(feature, "/mashr_meta_results.RData"))
    ## Estimate mixture proportions
    print(get_estimated_pi(m))
    plot_mixture_prop(m, feature)
    ## Save all lfsr
    data.frame(m$result$lfsr) %>%
        rename_all(list(~stringr::str_replace_all(., '.mean', ''))) %>%
        tibble::rownames_to_column("Effect") %>%
        data.table::fwrite(paste0(feature,"/lfsr_feature_4tissues.txt.gz"),
                           sep='\t')
    ## Save all effect size
    data.frame(m$result$PosteriorMean) %>%
        rename_all(list(~stringr::str_replace_all(., '.mean', ''))) %>%
        tibble::rownames_to_column("Effect") %>%
        data.table::fwrite(paste0(feature,"/posterior_mean_feature_4tissues.txt.gz"),
                           sep='\t')
}

## Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
args <- parser$parse_args()

## Run mashr for specific feature
run_mashr(args$feature)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
