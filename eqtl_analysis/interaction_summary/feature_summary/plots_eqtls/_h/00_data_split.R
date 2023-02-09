#### Script to subset expression
suppressMessages({
    library(argparse)
    library(magrittr)
})

#### MAIN analysis
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-t", "--tissue", type="character", default="Caudate",
                    help="Tissue to be analyzed [default: %default]")
args   <- parser$parse_args()

                                        # Set variables for filtering data
feature <- "Gene"; tissue  <- args$tissue
if(tissue == "Dentate Gyrus"){
    new_tissue <- "dentateGyrus"
} else {
    new_tissue <- gsub(" ", "_", tolower(tissue))
}

                                        # Generate output directory
outdir     <- paste(new_tissue, "out", sep="/")
dir.create(outdir, recursive=TRUE)

                                        # Extract interacting eQTL
infile  <- "../../_m/BrainSeq_ancestry_4features_4regions.txt.gz"
eqtls   <- read.table(gzfile(infile), comment.char="", header=TRUE, sep="\t") %>%
    dplyr::filter(Feature == feature, Tissue == tissue) %>%
    dplyr::arrange(gene_id, lfsr) %>% dplyr::group_by(gene_id) %>%
    dplyr::slice(1) %>% dplyr::arrange(lfsr) %>% as.data.frame
if(new_tissue == "dentateGyrus"){
    outfile <- paste0(outdir, "/eqtls.tsv")
    eqtls   <- eqtls
} else {
    outfile <- paste0(outdir, "/eqtls_top25.tsv")
    eqtls   <- eqtls %>% head(25)
}
eqtls %>% data.table::fwrite(outfile, sep="\t")

                                        # Load residuals
infile  <- here::here("eqtl_analysis/all_samples", new_tissue,
                      "covariates/residualized_expression/_m",
                      "genes_residualized_expression.csv")
dat     <- read.table(gzfile(infile), comment.char="", header=TRUE, sep=",") %>%
    dplyr::filter(gene_id %in% eqtls$gene_id)

                                        # Save data
write.table(dat, paste0(outdir,"/p"), col.names=TRUE,
            row.names=FALSE, quote=FALSE, sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
