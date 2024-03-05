## This script generates metric summary
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(stringr)
    library(argparse)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

read_metrics <- function(fname){
    feat_lt <- list("genes"="Gene", "transcripts"="Transcript",
                    "exons"="Exon", "junctions"="Junction")
    regions <- list("caudate"="Caudate", "dlpfc"="DLPFC",
                    "hippocampus"="Hippocampus",
                    "dentateGyrus"="Dentate Gyrus")
    ftype   <- str_split(str_extract(fname, "\\w+/r2.csv"),
                         "/", simplify=TRUE)[1]
    tissue  <- str_split(fname, "/", simplify=TRUE)[5]
    return(data.table::fread(fname) %>%
           mutate(type   = feat_lt[[ftype]],
                  region = regions[[tissue]]))
}

merge_data <- function(feature){
    fnames <- list.files("../../../_m", "r2.csv",
                         recursive=TRUE, full.names=TRUE)
    fnames <- fnames[grep(feature, fnames)]
    return(purrr::map(fnames, read_metrics) %>% bind_rows)
}

get_degs <- function(feature){
    feat_lt <- list("genes"="Gene", "transcripts"="Transcript",
                    "exons"="Exon", "junctions"="Junction")
    fname <- here::here("differential_analysis/tissue_comparison/summary_table/_m",
                        "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    return(data.table::fread(fname) %>%
           filter(Type == feat_lt[[feature]]) %>%
           select(Tissue, Effect, lfsr))
}

get_eGene <- function(feature){
    feat_lt <- list("genes"="Gene", "transcripts"="Transcript",
                    "exons"="Exon", "junctions"="Junction")
    fname <- here::here("eqtl_analysis/tissue_comparison/feature_summary/_m",
                        "BrainSeq_ancestry_4features_4regions.txt.gz")
    return(data.table::fread(fname) %>%
           filter(Feature == feat_lt[[feature]]) %>%
           select(Tissue, gene_id) %>% distinct %>%
           mutate(eGene = "eGene"))
}

plot_errplot <- function(dt, feature){
    dx     <- inner_join(dt, get_degs(feature), by=c("gene"="Effect",
                                                     "region"="Tissue")) %>%
        left_join(get_eGene(feature),
                  by=c("gene"="gene_id", "region"="Tissue")) %>%
        mutate(eGene = tidyr::replace_na(eGene, "Non-eGene"))
    errplt <-  dx %>% group_by(region, type) %>%
        mutate(popDE = ifelse(lfsr <= 0.05,
                              "Ancestry-Associated\nDEG (lfsr < 0.05)",
                              "Non-significant\nAncestry-Association Gene")) %>%
        ggerrorplot(x="region", y="r2_residual", color="popDE",
                    facet.by="eGene", palette="npg", ncol=1,
                    panel.labs.font=list(face='bold'),
                    ylab="Mean R2 (Test Score)", xlab="",
                    desc_stat="mean_ci", error.plot="pointrange",
                    ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=16) +
        rotate_x_text(45)
    save_plot(errplt, paste0("r2_correlation_", feature), 3.75, 8)
}

#### MAIN
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
args <- parser$parse_args()
                                        # Run analysis
feature <- args$feature
dt  <- merge_data(feature)
                                        # Summarize results
tb <- dt %>% filter(r2_residual > 0.01) %>%
    group_by(region, type) %>% count() %>%
    tidyr::pivot_wider(names_from=type, values_from=n)
print("Mean")
print(tb)
                                        # Plot metric
plot_errplot(dt, feature)
                                        # Save data
data.table::fwrite(dt,
                   file=paste0(feature,".prediction_metrics_summary.txt.gz"),
                   sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
