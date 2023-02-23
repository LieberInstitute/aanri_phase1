suppressMessages({
    library(here)
    library(dplyr)
    library(argparse)
})

#### MAIN analysis
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue to be analyzed [default: %default]")
parser$add_argument("-s", "--sge_id", type="integer",
                    help="Chunk to run analysis for")
args   <- parser$parse_args()
                                        # Run analysis
feature <- args$feature
tissue  <- args$tissue
sge_id  <- args$sge_id

                                        # Read feature annotation
indir  <- paste(tissue, feature, "out", sep="/")
pfile  <- paste0(indir, "/p", sge_id)
p      <- read.table(pfile, comment.char="", header=TRUE, sep="\t")

                                        # Select gene meta 
chr   <- gsub("chr","",p[,1])
wind  <- 200000
start <- ifelse(p[,2] - wind < 0, 0, p[,2] - wind)
end   <- p[,3] + wind
genes <- rownames(p)

outdir <- paste(tissue, feature, "ancestry", sep="/")
dir.create(outdir, showWarnings=FALSE)

for(j in 1:length(genes)){
    cat(j,"\n")
                                        # Get gene annotation
    chr1    <- chr[j]
    p1      <- start[j]
    p2      <- end[j]
    gname   <- genes[j]
if(feature == "junctions"){
    ## Clean up junction naming
    gname   <- gsub("\\(-\\)", "_minus",
                    gsub("\\(\\+\\)", "_plus", gname))
    gname   <- janitor::make_clean_names(gname)
    }
                                        # Run perl script
    outfile <- paste0(outdir, "/", gname)
    command <- paste("perl ../_h/local_ancestry_gene.pl",
                     "--chr", chr1, "--start", p1,
                     "--end", p2, "--out", outfile)
    system(command)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
