#### Script to subset expression
suppressMessages({
    library(argparse)
})

#### MAIN analysis
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue to be analyzed [default: %default]")
parser$add_argument("-c", "--chunks", type="integer", default=200,
                    help="Number of chunks to subset data [default: %default]")
args   <- parser$parse_args()
                                        # Run analysis
feature <- args$feature
tissue  <- args$tissue

infile  <- here::here("eqtl_analysis/",tissue,feature,
                      "normalize_expression/_m/",
                      paste0(feature, ".expression.bed.gz"))
dat     <- read.table(gzfile(infile), comment.char="", header=T)

                                        # output dir
outdir <- paste(tissue, feature,"out", sep="/")
dir.create(outdir, recursive=TRUE)
                                        # split dat
n  <- args$chunks
ps <- cut(seq(1,nrow(dat)),breaks=n,labels=FALSE)
for(i in 1:n){
    cat(i,"\n")
    idx <- which(ps == i,arr.ind=TRUE)
    outfile <- paste0(outdir,"/p",i)
    write.table(dat[idx,], outfile, col.names=TRUE,
                row.names=FALSE, quote=FALSE, sep="\t")
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
