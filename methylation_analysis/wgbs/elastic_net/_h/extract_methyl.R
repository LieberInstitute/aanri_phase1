#!/usr/bin/env Rscript
## @author = "Kynon JM Benjamin"

suppressMessages({
    library(here)
    library(dplyr)
    library(bsseq)
    library(argparse)
})

get_data <- function(chrom, tissue){
    flist <- list("hippocampus"=paste0("Hippocampus_chr",chrom,"BSobj_Genotypes.rda"),
                  "dlpfc"=paste0("dlpfc_chr",chrom,"BSobj_GenotypesDosage.rda"),
                  "caudate"=paste0("Caudate_chr",chrom,"BSobj_Genotypes.rda"))
                                        # Files
    ancestry <- here("input/ancestry_structure",
                     "structure.out_ancestry_proportion_raceDemo_compare")
    cpg_file <- here("input/methylation/_m/", flist[[tissue]])
    load(cpg_file, verbose=TRUE)
                                        # Phenotypes
    pheno <- colData(BSobj) %>% as.data.frame %>%
        inner_join(data.table::fread(ancestry),
                   by=c("BrNum"="id", "race"="group")) %>%
        filter(primarydx == "Control", race == "AA", agedeath > 17,
               Region == "Hippocampus")
                                        # Sample annotation
    brnum            <- pheno$BrNum
    rownames(pheno)  <- brnum
                                        # Coverage
    methyl <- getMeth(BSobj, type = "smooth") %>% as.data.frame %>%
        select(any_of(pheno$BrNum)) %>%
        mutate(cpg_id = paste0("cpg_",1:dim(BSobj)[1], "_chr",chrom)) %>%
        tibble::column_to_rownames("cpg_id")
                                        # CpG annotation
    cpg_annot <- granges(BSobj) %>% as.data.frame %>%
        mutate(cpg_id = paste0("cpg_",1:dim(BSobj)[1], "_chr",chrom)) %>%
        tibble::column_to_rownames("cpg_id")
    return(list("dna_m"=methyl, "annot"=cpg_annot))
}

## MAIN
parser  <- ArgumentParser()
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Brain region to extract methylation [default: %default]")
args    <- parser$parse_args()
feature <- args$feature

outdir <- "methylation"; dir.create(outdir)
for(chrom in c(1:22, "X")){
    lt <- get_data(chrom)
    fname1 <- paste0(outdir,"/methylation_wgbs_chr",chrom,".tsv")
    fname2 <- paste0(outdir,"/annotation_CpGs_wgbs_chr",chrom,".tsv")
    lt[["dna_m"]] %>% tibble::rownames_to_column("cpg_id") %>%
        data.table::fwrite(fname1, sep="\t")
    lt[["annot"]] %>% tibble::rownames_to_column("cpg_id") %>%
        data.table::fwrite(fname2, sep="\t")
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
