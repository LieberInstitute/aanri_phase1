## Calculated residualized expression values with covariates in eQTL analysis

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
})

                                        # Create parser object
option_list <- list(
    make_option(c("-f", "--feature"), type="character", default="genes",
                help="Feature to generate residualized expression [default: %default]",
                metavar="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

get_residualized <- function(feature){
    expr_file = paste0("../../../normalize_expression/_m/", feature,
                       ".expression.bed.gz")
    out_file = paste(feature, "residualized_expression.csv", sep="_")
    cov_file = paste0("../../_m/",feature,".combined_covariates.txt")
    expr = data.table::fread(expr_file) %>%
        tibble::column_to_rownames("gene_id") %>%
        select(starts_with("Br"))
    covs = data.table::fread(cov_file) %>%
        tibble::column_to_rownames("ID")
    mod = limma::lmFit(expr, t(rbind(intercept=1, covs[-1, ])))
    residuals(mod, expr) %>% as.data.frame %>%
        tibble::rownames_to_column("gene_id") %>%
        data.table::fwrite(out_file, sep=',')
}

if(length(opt$feature) == 0){
    print("Missing feature!")
} else {
    get_residualized(opt$feature)
}

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
