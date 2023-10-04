## This is for DMR annotation
suppressPackageStartupMessages({
    library(dplyr)
    library(argparse)
    library(annotatr)
    library(plyranges)
})

#### MAIN ####
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-t", "--tissue", type="character", default="dlpfc",
                    help="Tissue to be analyzed [default: %default]")
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
args   <- parser$parse_args()
tissue <- args$tissue; feature <- args$feature

                                        # Make directory
outdir <- file.path(tissue)
if(!dir.exists(outdir)){ dir.create(outdir, recursive=TRUE) }

                                        # Load all VMRs
fn  <- paste0("../../_m/aaOnly/vmr_", tissue, ".bed")
vmr <- data.table::fread(fn) |>
    rename("V1"="seqnames", "V2"="start", "V3"="end") |>
    mutate(vmr_id=paste(seqnames, start, end, sep="_")) |>
    as_granges()

                                        # Annotate with annotatr
                                        # Build the annotations
                                        # (a single GRanges object)
annotations <- build_annotations(genome='hg38',
                                 annotations=c('hg38_basicgenes'))

                                        # Intersect the regions we
                                        # read in with the annotations
vmr_annot   <- annotate_regions(
    regions=vmr,
    annotations=annotations,
    ignore.strand=TRUE,
    quiet=FALSE
)

                                        # A GRanges object is returned
outfile    <- paste0(outdir, "/vmr_annotation.tsv")
vmr_annot |> as.data.frame() |>
    data.table::fwrite(outfile, sep='\t')

vmr_annsum <- summarize_annotations(
    annotated_regions=vmr_annot,
    quiet=TRUE
)
print(vmr_annsum)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
