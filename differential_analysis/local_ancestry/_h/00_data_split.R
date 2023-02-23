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
config  <- list("genes"="gene", "transcripts"="tx",
                "exons"="exon", "junctions"="jxn")
                                        # Load and filter feature annotation
load(here::here("differential_analysis", tissue, "_m",
                feature,"voomSVA.RData"))
infile  <- here::here("input/text_files_counts/_m", tissue,
                      paste0(config[feature],"_annotation.tsv"))
dat     <- data.table::fread(infile) |>
    tibble::column_to_rownames("names")
dat     <- dat[rownames(v$genes), ]
rm(v)
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
                row.names=TRUE, quote=FALSE, sep="\t")
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
