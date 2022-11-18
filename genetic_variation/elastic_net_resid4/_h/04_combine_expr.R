## This script combines files
suppressPackageStartupMessages({
    library(argparse)
    library(tidyverse)
})

read_data <- function(fname){
    return(read.csv(fname, row.names=1))
}

#### MAIN
                                        # Create parser object
parser  <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to combine results [default: %default]")
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue to combine results [default: %default]")
args    <- parser$parse_args()
feature <- args$feature; tissue <- args$tissue

                                        # Residualized Expression
outfile1 <- paste0(feature, ".predicted_expression.",
                   tissue, ".txt.gz")
flist1 <- list.files(paste(tissue, feature, "out", sep="/"),
                     pattern="res_pred.csv", full.names=TRUE)
map(flist1, read_data) %>% bind_cols %>% t %>% as.data.frame %>%
    rownames_to_column("gene_id") %>% na.omit %>%
    data.table::fwrite(outfile1, sep='\t', compress="gzip")

                                        # Normalized Expression
outfile2 <- paste0(feature, ".predicted_expression_norm.",
                   tissue, ".txt.gz")
flist2 <- list.files(paste(tissue, feature, "out", sep="/"),
                     pattern="norm_pred.csv", full.names=TRUE)
map(flist2, read_data) %>% bind_cols %>% t %>% as.data.frame %>%
    rownames_to_column("gene_id") %>% na.omit %>%
    data.table::fwrite(outfile2, sep='\t', compress="gzip")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
