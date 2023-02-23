## This script combines files
suppressPackageStartupMessages({
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

                                        # Residualized Expression
outfile <- paste0(feature,"/diffExpr_local.tsv")
flist   <- list.files(feature, pattern="\\.ancestry", full.names=TRUE)
purrr::map(flist, read_data) |> dplyr::bind_rows() |>
    data.table::fwrite(outfile, sep='\t')

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
