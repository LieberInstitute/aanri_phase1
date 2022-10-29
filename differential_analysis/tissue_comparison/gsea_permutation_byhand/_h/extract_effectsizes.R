## This script randomly assigns variable of interest and performs
## limma-voom and mash modeling to get effect sizes. It will
## require that entrez or ensembl ids are provided in the annotation
## file.

library(mashr)
library(optparse)

#### Functions
prep_data <- function(perm_num){
    permutation_file = paste0("_model_bhat_N_shat_perm",perm_num,".txt")
    files_lt <- list.files("temp_files", pattern=paste0("*",permutation_file))
    bhat_list <- list(); shat_list <- list()
    for(ii in seq_along(files_lt)){
        fn    <- files_lt[ii]
        label <- regmatches(fn, regexpr("^[^_]+(?=_)", fn, perl=TRUE))
        df    <- read.table(paste0("temp_files/",fn), sep='\t', header=TRUE)
        bhat  <- data.frame(df$feature_ids, df$bhat);
        shat  <- data.frame(df$feature_ids, df$shat)
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

run_mashr <- function(perm_num){
    ## Load prepared data
    datalist <- prep_data(perm_num)
    bhat <- as.matrix(datalist[["bhat"]])
    shat <- as.matrix(datalist[["shat"]])
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
    ## Save all effect size
    df <- data.frame(m$result$PosteriorMean)
    write.table(df, sep="\t", row.names=TRUE, quote=FALSE,
                file=paste0("temp_files/mash_effectsize_perm",perm_num,".txt"))
}

#### MAIN ####
                                        # parser
option_list <- list(
    make_option(c("-p", "--perm_num"), type="integer", default=1,
                help="permutation number [default=%default]")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)
                                        # Run mash
run_mashr(opt$perm_num)
