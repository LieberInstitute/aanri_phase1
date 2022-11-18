library(argparse)

#### MAIN analysis
                                        # Create parser object
parser  <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue to be analyzed [default: %default]")
args    <- parser$parse_args()
feature <- args$feature
tissue  <- args$tissue

                                        # results of residualized expression
files <- list.files(paste(tissue, feature, "out", sep="/"),
                    "res.csv", full.names=TRUE)
dat   <- dplyr::bind_rows(purrr::map(files, read.csv))
m2    <- apply(dat[,-(1:2)],1,function(x){mean(x^2,na.rm=T)})

                                        # combine results and save
res     <- data.frame(gene=dat$genes,r2_residual=m2)
outfile <- paste(tissue, feature, "r2.csv", sep="/")
write.csv(res, outfile, row.names=FALSE)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
