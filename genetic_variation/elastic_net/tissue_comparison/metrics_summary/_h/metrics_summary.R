## This script generates metric summary
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
})

read_metrics <- function(fname){
    feat_lt <- list("genes"="Gene", "transcripts"="Transcript",
                    "exons"="Exon", "junctions"="Junction")
    regions <- list("caudate"="Caudate", "dlpfc"="DLPFC",
                    "hippocampus"="Hippocampus",
                    "dentateGyrus"="Dentate_Gyrus")
    ftype   <- str_split(str_extract(fname, "\\w+.pred"),
                         "\\.", simplify=TRUE)[1]
    tissue  <- str_split(fname, "\\.", simplify=TRUE)[9]
    return(data.table::fread(fname) %>% group_by(feature) %>%
           summarise(across(n_features:test_score_evar,
                            list(median = median, mean = mean,
                                 std = sd))) %>%
           as.data.frame %>%
           mutate(type   = feat_lt[[ftype]],
                  region = regions[[tissue]]))
}

merge_data <- function(){
    fnames <- list.files("../../../_m", "*prediction_metrics*",
                         full.names=TRUE)
    return(purrr::map(fnames, read_metrics) %>% bind_rows)
}

#### MAIN
dt  <- merge_data()
                                        # Summarize results
tb1 <- dt %>% filter(test_score_r2_mean > 0.01) %>%
    group_by(region, type) %>% count() %>%
    tidyr::pivot_wider(names_from=type, values_from=n) %>%
    select(Gene, Transcript, Exon, Junction)
print("Mean")
print(tb1)
tb2 <- dt %>% filter(test_score_r2_median > 0.01) %>%
    group_by(region, type) %>% count() %>%
    tidyr::pivot_wider(names_from=type, values_from=n) %>%
    select(Gene, Transcript, Exon, Junction)
print("Median")
print(tb2)
                                        # Save data
data.table::fwrite(dt, file="prediction_metrics_summary.txt.gz",
                   sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
