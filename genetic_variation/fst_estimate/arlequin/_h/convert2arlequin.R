## This script is used to convert PLINK to ARLEQUIN format.
library(dplyr)
library(radiator)

generate_strata <- function(){
    strata_file <- "../../_m/within_clusters.txt"
    strata <- data.table::fread(strata_file) %>% select(FID, group)
    colnames(strata) <- c("INDIVIDUALS", "STRATA")
    strata %>% data.table::fwrite("strata.tsv", sep="\t")
}

radiator::summary_strata("strata.tsv")

radiator::detect_genomic_format(data="../../_m/caudate/ENSG00000111664.10.bed")

df <- radiator::genomic_converter(data="../../_m/caudate/ENSG00000111664.10.bed",
                                  strata = "strata.tsv",
                                  output = c("genepop", "arlequin"),
                                  filename = "ENSG00000111664.10",
                                  parallel.core = parallel::detectCores() - 1,
                                  verbose = TRUE)
genepop <- df$genepop
arlequin <- df$arlequin
