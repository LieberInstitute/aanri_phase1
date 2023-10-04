## This script combines files
suppressPackageStartupMessages({
    library(dplyr)
    library(argparse)
})

read_data <- function(fname){
    return(data.table::fread(fname))
}

#### MAIN
                                        # Create parser object
parser  <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to combine results [default: %default]")
args    <- parser$parse_args()
feature <- args$feature

## Generate files
##                                         # Main model
## outfile1 <- paste0(feature,"/diffExpr_model_main.tsv")
## flist1   <- list.files(feature, pattern="\\.main", full.names=TRUE)
## purrr::map(flist1, read_data) |> bind_rows() |>
##     data.table::fwrite(outfile1, sep='\t')

##                                         # VMR model
## outfile2 <- paste0(feature,"/diffExpr_model_vmr.tsv")
## flist2   <- list.files(feature, pattern="\\.vmr", full.names=TRUE)
## purrr::map(flist2, read_data) |> bind_rows() |>
##     data.table::fwrite(outfile2, sep='\t')
                                        # Partial r2
outfile <- paste0(feature,"/partial_r2.tsv")
                                        # Main model
flist1  <- list.files(feature, pattern="\\.general.r2", full.names=TRUE)
gen_df  <- purrr::map(flist1, read_data) |> bind_rows()
                                        # VMR model
flist2  <- list.files(feature, pattern="\\.ancestry.r2", full.names=TRUE)
vmr_df <- purrr::map(flist2, read_data) |> bind_rows() 
                                        # Save data
inner_join(gen_df, vmr_df, by="feature_id", suffix=c(".gen", ".vmr")) |>
    mutate(partial_r2=(partial_r2.gen - partial_r2.vmr) / partial_r2.gen) |>
    arrange(desc(partial_r2)) |>
    select(feature_id, partial_r2, partial_r2.gen, partial_r2.vmr) |>
    data.table::fwrite(outfile, sep='\t')

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
