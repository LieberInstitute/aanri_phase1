## This script combines files
suppressPackageStartupMessages({
    library(argparse)
    library(tidyverse)
})

read_data <- function(fname){
    return(data.table::fread(fname) %>%
           column_to_rownames("fid"))
}

read_metrics <- function(fname){
    return(data.table::fread(fname))
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
                                        # Expression
outfile1 <- paste0(feature, ".predicted_expression.",
                   tissue, ".txt.gz")
flist1 <- list.files(paste(tissue, feature, sep="/"),
                     pattern="^pred*", full.names=TRUE)
map(flist1, read_data) %>% bind_cols %>% t %>%
    as.data.frame %>% rownames_to_column("gene_id") %>%
    data.table::fwrite(outfile1, sep='\t', compress="gzip")
                                        # Metrics
outfile2 <- paste0(feature, ".prediction_metrics.",
                   tissue, ".txt.gz")
flist2 <- list.files(paste(tissue, feature, sep="/"),
                     pattern="^enet*", full.names=TRUE)
map(flist2, read_metrics) %>% bind_rows %>%
    data.table::fwrite(outfile2, sep='\t', compress="gzip")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
