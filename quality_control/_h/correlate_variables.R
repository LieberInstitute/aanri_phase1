## Select and filter data for MDD analysis
suppressPackageStartupMessages({
    library(dplyr)
    library(argparse)
    library(SummarizedExperiment)
})

                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-r", "--region", type="character", default="caudate",
                    help="Brain region to extract [default: %default]")
args <- parser$parse_args()

                                        # Function from jaffelab github
merge_rse_metrics <- function(rse) {
    stopifnot(is(rse, 'RangedSummarizedExperiment'))
    rse$numMapped = sapply(rse$numMapped, sum)
    rse$numReads = sapply(rse$numReads, sum)
    rse$numUnmapped = sapply(rse$numUnmapped, sum)
    rse$mitoMapped = sapply(rse$mitoMapped, sum)
    rse$totalMapped = sapply(rse$totalMapped, sum)
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
    return(rse)
}

get_mds <- function(){
    baseloc  <- "/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/input"
    mds_file <- paste0(baseloc,"/genotypes/old_mds/_m/LIBD_Brain_TopMed.mds")
    mds = data.table::fread(mds_file) %>%
        rename_at(.vars = vars(starts_with("C")),
                  function(x){sub("C", "snpPC", x)}) %>%
        mutate_if(is.character, as.factor)
    return(mds)
}
memMDS <- memoise::memoise(get_mds)

prep_data <- function(region){
    baseloc  <- "/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/input"
    ancestry <- paste0(baseloc,"/ancestry_structure/structure.out_ancestry",
                       "_proportion_raceDemo_compare")
    counts_lt = list("caudate"=paste0(baseloc,"/counts/_m/caudate_brainseq_",
                                      "phase3_hg38_rseGene_merged_n464.rda"), 
                     "dentateGyrus"=paste0(baseloc,"/counts/_m/astellas_dg_",
                                           "hg38_rseGene_n263.rda"),
                     "dlpfc"=paste0(baseloc,"/counts/_m/dlpfc_ribozero_brainseq",
                                    "_phase2_hg38_rseGene_merged_n453.rda"),
                     "hippocampus"=paste0(baseloc,"/counts/_m/hippo_brainseq_",
                                          "phase2_hg38_rseGene_merged_n447.rda"))
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

remove_variables <- function(pheno_df){
    if(length(check_dup(pheno_df)) != 0){
        pheno_df <- pheno_df %>% select(-check_dup(pheno_df))
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
#gyrus   <- remove_variables(gyrus) ## by hand, drop only concordMapRate
dlpfc   <- remove_variables(dlpfc)
hippo   <- remove_variables(hippo)
print(cor.test(gyrus$concordMapRate, gyrus$overallMapRate))

## Limitation is gyrus available data

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
