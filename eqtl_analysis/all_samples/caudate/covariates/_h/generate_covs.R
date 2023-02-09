## Prepare covariates for FastQTL / tensorQTL analysis
suppressPackageStartupMessages({
    library(argparse)
    library(tidyverse)
})
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to extract [default: %default]")
parser$add_argument("-r", "--region", type="character", default="caudate",
                    help="Brain region to extract [default: %default]")
args <- parser$parse_args()

                                        # Load data
sample_df = data.table::fread("../../_m/sample_id_to_brnum.tsv")
covs_file <- here::here("input/covariates/all_samples", args$region,
                        "_m/covariates.rda")
load(covs_file, verbose=TRUE)
                                        # Configure
covs_lt = list("genes"=covsGene0, "transcripts"=covsTx0,
               "exons"=covsExon0, "junctions"=covsJxn0)
covs0 = covs_lt[[args$feature]]
                                        # Extract covariates
covs = covs0 %>% as.data.frame %>% select(sample_df$RNum) %>%
    t %>% as.data.frame %>% rownames_to_column("RNum") %>%
    inner_join(sample_df, by=c("RNum")) %>% select(-"RNum") %>%
    rename("ID"="BrNum") %>% column_to_rownames("ID") %>% t %>%
    as.data.frame %>% rownames_to_column("ID")
                                        # Save file
covs %>% data.table::fwrite(paste0(args$feature, ".combined_covariates.txt"),
                            sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
