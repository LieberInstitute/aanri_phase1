## This script prepares and runs mash modeling.

library(mashr)
library(optparse)

#### Functions
prep_data <- function(feature){
    de_file <- paste0("diffExpr_EAvsAA_full.txt")
    files_lt  <- list.files(c(paste0("../../caudate/_m/", feature),
                              paste0("../../dentateGyrus/_m/", feature),
                              paste0("../../dlpfc/_m/", feature),
                              paste0("../../hippocampus/_m/", feature)),
                            pattern=de_file, full.names=TRUE)
    bhat_list <- list(); shat_list <- list()
    for(ii in seq_along(files_lt)){
        fn    <- files_lt[ii]
        label <- basename(regmatches(fn, regexpr("^[^_]+(?=_)", fn, perl=TRUE)))
        df    <- read.table(fn, sep='\t', header=TRUE) |>
            dplyr::rename(feature_id=X, bhat=logFC, shat=SE)
        bhat  <- data.frame(df$feature_id, df$bhat);
        shat  <- data.frame(df$feature_id, df$shat)
        colnames(bhat) <- c("feature_id", label);
        colnames(shat) <- c("feature_id", label)
        bhat_list[[ii]] <- bhat; shat_list[[ii]] <- shat
    }
    bhat_dt <- Reduce(function(...) merge(..., by="feature_id"), bhat_list)
    rownames(bhat_dt) <- bhat_dt$feature_id
    shat_dt <- Reduce(function(...) merge(..., by="feature_id"), shat_list)
    rownames(shat_dt) <- shat_dt$feature_id
    return(list("bhat"=bhat_dt[, which(colnames(bhat_dt) != "feature_id")],
                "shat"=shat_dt[, which(colnames(shat_dt) != "feature_id")]))
}

run_mashr <- function(feature){
    ## Load prepared data
    datalist <- prep_data(feature)
    bhat     <- as.matrix(datalist[["bhat"]])
    shat     <- as.matrix(datalist[["shat"]])
    ## Prepare for mashr
    data.init = mash_set_data(bhat, shat)
    Vhat = estimate_null_correlation_simple(data.init)
    data = mash_update_data(data.init, V=Vhat)
    rm(data.init)
    ## Calculate data driven covariances
    U.c = cov_canonical(data)
    m.1by1 = mash_1by1(data)
    strong = get_significant_results(m.1by1)
    rm(m.1by1)
    U.pca = cov_pca(data, dim(datalist[["bhat"]])[2], subset=strong)
    U.ed = cov_ed(data, U.pca, subset=strong)
    rm(strong, U.pca)
    ## Fit mash model and compute posterior summaries
    m = mash(data, Ulist = c(U.ed,U.c), algorithm.version="Rcpp")
    rm(U.ed, U.c)
    ## Save all lfsr
    save(m, data, file=paste0(feature, "/mashr_meta_results.RData"))
    data.frame(m$result$lfsr) |>
        tibble::rownames_to_column("feature_id") |>
        write.table(sep="\t", row.names=FALSE, quote=FALSE,
                    file=paste0(feature, "/mash_lfsr.txt"))
    ## Save all effect size
    data.frame(m$result$PosteriorMean) |>
        tibble::rownames_to_column("feature_id") |>
        write.table(sep="\t", row.names=FALSE, quote=FALSE,
                    file=paste0(feature, "/mash_effectsize.txt"))
}

#### MAIN ####
                                        # Parser
option_list <- list(
    make_option(c("--feature"), type="character", default="genes",
                help="feature to analyze [default=%default]",
                metavar="character")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)
                                        # Run mash
dir.create(opt$feature)
run_mashr(opt$feature)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
