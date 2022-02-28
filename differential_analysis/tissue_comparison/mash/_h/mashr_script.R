suppressPackageStartupMessages({
    library(mashr)
    library(dplyr)
    library(ggpubr)
    library(argparse)
})

save_img <- function(image, fn, w, h){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_pvals <- function(feature){
    fn = paste0(feature, "/pvalue_de_4tissues.tsv")
    return(data.table::fread(fn) %>% mutate(effect=Feature))
}

get_bhat <- function(feature){
    fn = paste0(feature, "/bhat_de_4tissues.tsv")
    return(data.table::fread(fn) %>% mutate(effect=Feature) %>%
           distinct(effect, .keep_all=TRUE))
}

get_shat <- function(feature){
    fn = paste0(feature, "/shat_de_4tissues.tsv")
    return(data.table::fread(fn) %>% mutate(effect=Feature) %>%
           distinct(effect, .keep_all=TRUE))
}

save_results <- function(pval, m, feature){
    fn = paste0(feature,"/brainseq_ancestry_4tissues_mashr.tsv")
    sig_de <- get_significant_results(m)
    df = get_n_significant_conditions(m, thresh = 0.05, conditions = NULL,
                                      sig_fn = get_lfsr)
    cc = get_n_significant_conditions(m, thresh = 0.05, conditions = 1,
                                      sig_fn = get_lfsr)
    gg = get_n_significant_conditions(m, thresh = 0.05, conditions = 2,
                                      sig_fn = get_lfsr)
    dd = get_n_significant_conditions(m, thresh = 0.05, conditions = 3,
                                      sig_fn = get_lfsr)
    hh = get_n_significant_conditions(m, thresh = 0.05, conditions = 4,
                                      sig_fn = get_lfsr)
    dt = data.frame(N_Regions_Shared=df, Caudate=cc, `Dentate Gyrus`=gg,
                    DLPFC=dd, Hippocampus=hh) %>%
        tibble::rownames_to_column("effect")
    print("Sharing across brain regions:")
    print(table(dt$N_Regions_Shared))
    ## Total number of tissue specific eQTL
    print(paste("There are",sum(df == 1), feature, "with specific expression!"))
    dt %>% rename("Dentate Gyrus"="Dentate.Gyrus", "Feature"="effect") %>%
        data.table::fwrite(fn, sep='\t')
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
    bhat <- get_bhat(feature) %>% tibble::column_to_rownames("effect") %>%
        select(-Feature) %>% as.matrix
    shat <- get_shat(feature) %>% tibble::column_to_rownames("effect") %>%
        select(-Feature) %>% as.matrix
    pval <- get_pvals(feature)
    ## Prepare for mashr
    data.init = mash_set_data(bhat, shat)
    Vhat = estimate_null_correlation_simple(data.init)
    rm(data.init)
    ## Initial mash
    data = mash_set_data(bhat,shat, V=Vhat)
    ## Calculate data driven covariances
    U.pca = cov_pca(data,4)
    U.ed = cov_ed(data, U.pca)
    U.c = cov_canonical(data)
    ## Fit mash model and compute posterior summaries
    m = mash(data, Ulist = c(U.ed,U.c))
    ## Examine pairwise sharing
    print(get_pairwise_sharing(m))
    ## Save significant results
    print(get_significant_results(m) %>% length)
    save_results(pval, m, feature)
    save(m, file=paste0(feature, "/mashr_meta_results.RData"))
    ## Estimate mixture proportions
    print(get_estimated_pi(m))
    plot_mixture_prop(m, feature)
    ## Save all lfsr
    data.frame(m$result$lfsr) %>%
        tibble::rownames_to_column("Effect") %>%
        data.table::fwrite(paste0(feature,"/lfsr_feature_4tissues.txt.gz"),
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
