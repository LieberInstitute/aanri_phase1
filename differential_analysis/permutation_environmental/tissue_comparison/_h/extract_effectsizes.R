## This script prepares and runs mash modeling.

library(mashr)
library(optparse)

#### Functions
prep_data <- function(feature, perm){
    perm_num  <- stringr::str_pad(perm, 2, pad="0")
    perm_file <- paste0("de_model_ancestry_perm_",perm_num,".txt")
    files_lt  <- list.files(c(paste0("../../caudate/_m/", feature),
                              paste0("../../dentateGyrus/_m/", feature),
                              paste0("../../dlpfc/_m/", feature),
                              paste0("../../hippocampus/_m/", feature)),
                            pattern=perm_file, full.names=TRUE)
    bhat_list <- list(); shat_list <- list()
    for(ii in seq_along(files_lt)){
        fn    <- files_lt[ii]
        label <- basename(regmatches(fn, regexpr("^[^_]+(?=_)", fn, perl=TRUE)))
        df    <- read.table(fn, sep='\t', header=TRUE)
        bhat  <- data.frame(df$X, df$logFC); shat  <- data.frame(df$X, df$SE)
        colnames(bhat) <- c("ids", label); colnames(shat) <- c("ids", label)
        bhat_list[[ii]] <- bhat; shat_list[[ii]] <- shat
    }
    bhat_dt <- Reduce(function(...) merge(..., by="ids"), bhat_list)
    rownames(bhat_dt) <- bhat_dt$ids
    shat_dt <- Reduce(function(...) merge(..., by="ids"), shat_list)
    rownames(shat_dt) <- shat_dt$ids
    return(list("bhat"=bhat_dt[, which(colnames(bhat_dt) != "ids")],
                "shat"=shat_dt[, which(colnames(shat_dt) != "ids")]))
}

run_mashr <- function(feature, perm){
    perm_num <- stringr::str_pad(perm, 2, pad="0")
    ## Load prepared data
    datalist <- prep_data(feature, perm)
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
    df <- tibble::rownames_to_column(data.frame(m$result$lfsr), "feature_id")
    write.table(df, sep="\t", row.names=FALSE, quote=FALSE,
                file=paste0(feature, "/mash_lfsr_perm_",perm_num,".txt"))
    ## Save all effect size
    df <- tibble::rownames_to_column(data.frame(m$result$PosteriorMean),
                                     "feature_id")
    write.table(df, sep="\t", row.names=FALSE, quote=FALSE,
                file=paste0(feature, "/mash_effectsize_perm_",perm_num,".txt"))
}

#### MAIN ####
                                        # Parser
option_list <- list(
    make_option(c("--feature"), type="character", default="genes",
                help="feature to analyze [default=%default]",
                metavar="character"),
    make_option(c("-p", "--perm_num"), type="integer", default=1,
                help="permutation number [default=%default]")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)
                                        # Run mash
dir.create(opt$feature)
run_mashr(opt$feature, opt$perm_num)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
