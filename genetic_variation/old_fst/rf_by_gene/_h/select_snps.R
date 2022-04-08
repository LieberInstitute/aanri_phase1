## Calculate the partial R2 using RaFFE SNPs and individually

suppressPackageStartupMessages({
    library("argparse")
    library("tidyverse")
})

get_ml_summary <- function(Tissue, geneid, model){
    fn = paste0("../../../snp_prediction/de_genes/", model,
                "/summary_10Folds_allTissues.tsv")
    ml_df = data.table::fread(fn) %>% mutate_at("fold", as.character) %>%
        filter(tissue == Tissue, feature==geneid) %>%
        select(fold, n_features, n_redundant, starts_with("test_score_r2")) %>%
        pivot_longer(-fold) %>% group_by(name) %>%
        summarise(Mean=mean(value), Median=median(value), Std=sd(value),
                  .groups = "keep") %>% as.data.frame %>%
        column_to_rownames("name")
    return(ml_df)
}

get_snps_ml <- function(tissue, geneid, model){
    fn = paste0("../../../snp_prediction/de_genes/", model,
                "/rank_features_allTissues.tsv")
    num <- get_ml_summary(tissue, geneid, model)["n_features", "Median"]
    snps = data.table::fread(fn) %>% filter(Tissue==tissue, Feature==geneid) %>%
        group_by(Fold, SNP) %>% select(-c(Tissue, Feature)) %>% slice(1) %>%
        arrange(Fold, desc(Rank))%>% as.data.frame %>%
        pivot_wider(names_from=Fold, values_from=Rank) %>%
        mutate(median_all = apply(., 1, median)) %>% arrange(median_all) %>%
        mutate(Rank = rank(median_all)) %>% filter(Rank <= num) %>% select(SNP)
    return(snps)
}

tissue_map <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus",
                "dentateGyrus"="Dentate Gyrus")[[tissue]])
}

save_snps <- function(tissue, geneid, model){
    gname = gsub("\\.", "_", geneid)
    snps = get_snps_ml(tissue_map(tissue), gname, model)
    unique_snps <- str_split(snps$SNP, "_", simplify=TRUE)[, 1] %>% unique
    data.frame("SNP"=unique_snps) %>%
        data.table::fwrite(paste0(tissue, "/", geneid, "_snps.txt"),
                           sep='\t', col.names=FALSE)
}

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-g", "--gname", type="character",
                    help="Name of gene to analyze [required]")
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue for brain analysis [default: %default]")
parser$add_argument("-m", "--model", type="character", default="rf",
                    help="ML model for directory [default: %default]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

if(length(args$gname) == 0){
    print("Missing gene to analyze!")
} else {
    save_snps(args$tissue, args$gname, args$model)
}
