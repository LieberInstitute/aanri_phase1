## Select and filter data for MDD analysis
suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(SummarizedExperiment)
})
merge_rse_metrics <- function(rse) {
                                        # Function from jaffelab github
    stopifnot(is(rse, 'RangedSummarizedExperiment'))
    rse$overallMapRate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$overallMapRate, rse$numReads)
    rse$mitoRate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$mitoRate, rse$numMapped)
    rse$rRNA_rate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$rRNA_rate, rse$numMapped)
    rse$totalAssignedGene = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$totalAssignedGene, rse$numMapped)
    rse$numMapped = sapply(rse$numMapped, sum)
    rse$numReads = sapply(rse$numReads, sum)
    rse$numUnmapped = sapply(rse$numUnmapped, sum)
    rse$mitoMapped = sapply(rse$mitoMapped, sum)
    rse$totalMapped = sapply(rse$totalMapped, sum)
    return(rse)
}

get_mds <- function(){
    mds_file <- here("input/genotypes/old_mds/_m/LIBD_Brain_TopMed.mds")
    mds = data.table::fread(mds_file) %>%
        rename_at(.vars = vars(starts_with("C")),
                  function(x){sub("C", "snpPC", x)}) %>%
        mutate_if(is.character, as.factor)
    return(mds)
}
memMDS <- memoise::memoise(get_mds)

prep_data <- function(region){
    ancestry <- here("input/ancestry_structure/",
                     "structure.out_ancestry_proportion_raceDemo_compare")
    counts_lt = list("caudate"=here("input/counts/_m/",
                                    paste0("caudate_brainseq_phase3_hg38",
                                           "_rseGene_merged_n464.rda")),
                     "dentateGyrus"=here("input/counts/_m/",
                                         "astellas_dg_hg38_rseGene_n263.rda"),
                     "dlpfc"=here("input/counts/_m/",
                                  paste0("dlpfc_ribozero_brainseq_phase2",
                                         "_hg38_rseGene_merged_n453.rda")),
                     "hippocampus"=here("input/counts/_m/",
                                        paste0("hippo_brainseq_phase2_hg38",
                                               "_rseGene_merged_n447.rda")))
    load(counts_lt[[region]], verbose=TRUE)
    rse_df = rse_gene
    keepIndex = which((rse_df$Dx %in% c("Control", "Schizo")) & 
                      rse_df$Age > 17 & rse_df$Race %in% c("AA", "CAUC"))
    rse_df = rse_df[, keepIndex]
    rse_df$Dx = factor(rse_df$Dx, levels = c("Control", "Schizo"))
    rse_df$Sex <- factor(rse_df$Sex)
    rse_df <- merge_rse_metrics(rse_df)
    colData(rse_df)$RIN = sapply(colData(rse_df)$RIN,"[",1)
    rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"), "[", 1)
    pheno = colData(rse_df) %>% as.data.frame %>% 
        select("RNum", "BrNum", "Dx", "Race", "Sex", "Age", "RIN", "mitoRate",
               "totalAssignedGene", "overallMapRate", "concordMapRate","rRNA_rate",
               any_of(starts_with(c("Batch", "Adapter", "phredGT", "percent"))),
               -any_of("AdapterContent")) %>%
        mutate_if(is.list, ~sapply(., sum)) %>%
        mutate_if(is.numeric, scales::rescale) %>%
        inner_join(data.table::fread(ancestry), by=c("BrNum"="id", "Race"="group")) %>%
        inner_join(memMDS() %>% select("FID", "snpPC1", "snpPC2", "snpPC3"),
                   by=c("BrNum"="FID")) %>% distinct(RNum, .keep_all = TRUE)
    return(pheno)
}

check_dup <- function(df){
    sample <- df %>% select_if(is.numeric)
    variables <- names(sample)
    return(cytominer::correlation_threshold(variables, sample, cutoff=0.90))
}

check_corr <- function(df){
    sample <- df %>% select_if(is.numeric)
    dt = sample %>% corrr::correlate() %>%
        corrr::stretch() %>% tidyr::drop_na() %>%
        filter(abs(r) > 0.90) %>%
        distinct(r, .keep_all=TRUE)
    varX <- distinct(dt, x)$x
    varX <- varX[-which(varX %in% intersect(varX, distinct(dt, y)$y))]
    vars <- unique(c(distinct(dt, x)$x, distinct(dt, y)$y))
    return(setdiff(vars, varX))
}

old_remove_variables <- function(pheno_df){
    if(length(check_dup(pheno_df)) != 0){
        pheno_df <- pheno_df %>% select(-check_dup(pheno_df))
    }
    return(pheno_df)
}

remove_variables <- function(pheno_df){
    if(length(check_corr(pheno_df)) != 0){
        pheno_df <- pheno_df %>% select(-check_corr(pheno_df))
    }
    return(pheno_df)
}

#### MAIN ####
                                        # Load phenotypes
caudate <- prep_data("caudate")
gyrus   <- prep_data("dentateGyrus")
dlpfc   <- prep_data("dlpfc")
hippo   <- prep_data("hippocampus")
                                        # Drop correlated
caudate <- remove_variables(caudate)
gyrus   <- remove_variables(gyrus)
dlpfc   <- remove_variables(dlpfc)
hippo   <- remove_variables(hippo)
                                        # Commone variables
vars <- intersect(colnames(caudate),
                  intersect(colnames(gyrus),
                            intersect(colnames(dlpfc), colnames(hippo))))
data.frame("Variables"=vars) %>%
    data.table::fwrite("shared_variables.tsv", sep='\t')
                                        # Save variables that
                                        # are not highly correlated
caudate %>% data.table::fwrite("caudate_expression_covs.tsv", sep='\t')
gyrus %>% data.table::fwrite("dentateGyrus_expression_covs.tsv", sep='\t')
dlpfc %>% data.table::fwrite("dlpfc_expression_covs.tsv", sep='\t')
hippo %>% data.table::fwrite("hippocampus_expression_covs.tsv", sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
