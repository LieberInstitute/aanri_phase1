#### This script calculates the absolute allele frequency differences
#### between 1000 Genome Project Afr and Eur references genome per
#### chromosome.

load_afr <- function(chrom){
    fn <- here::here("input/genotypes/gp1000/_m",
                     paste0("AFR_noIndels_chr", chrom, "_maf0.01.frq"))
    return(data.table::fread(fn))
}

load_eur <- function(chrom){
    fn <- here::here("input/genotypes/gp1000/_m",
                     paste0("EUR_noIndels_chr",chrom,"_maf0.005.frq"))
    return(data.table::fread(fn))
}

merge_data <- function(chrom){
    return(dplyr::inner_join(load_afr(chrom), load_eur(chrom),
                             by=c("CHR", "SNP"),
                             suffix=c("_AFR", "_EUR")))
}

cal_allele_diff <- function(chrom){
    return(merge_data(chrom) |>
           dplyr::mutate(AFD=ifelse(A1_AFR == A1_EUR,
                                    abs(MAF_AFR - MAF_EUR),
                                    abs(MAF_AFR - (1-MAF_EUR)))) |>
           dplyr::select(CHR, SNP, A1_AFR, A2_AFR, AFD))
}

#### MAIN
for(chrom in seq(1, 22, 1)){
    outfile <- paste0("allele_freq_difference.1000GP.chr", chrom, ".tsv")
    data.table::fwrite(cal_allele_diff(chrom), outfile, sep="\t")
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

