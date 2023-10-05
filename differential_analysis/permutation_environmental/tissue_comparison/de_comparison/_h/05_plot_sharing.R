## Generate plots comparing Black only and Black v. White
## results venn diagram showing overlap

suppressPackageStartupMessages({
    library(dplyr)
    library(ggvenn)
    library(ggplot2)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_comp_data <- function(feature, tissue){
    fn <- "BrainSeq_ancestry_DE_comparison.txt.gz"
    return(data.table::fread(fn) |>
           filter(region == tissue, feature_type == feature))
}
memCOMP <- memoise::memoise(load_comp_data)

venn_diagrams <- function(feature, tissue){
    new_feature <- tolower(feature)
    new_tissue  <- gsub(" ", "_", tolower(tissue))
    outfile     <- paste("venn_diagram", new_feature,
                         new_tissue, sep=".")
    x  <- list(
        `Black only`     = memCOMP(feature, tissue) |>
            filter(ancestry == 1) |> pull(feature_id),
        `Black v. White` = memCOMP(feature, tissue) |>
            filter(environ == 1) |> pull(feature_id)
    )
    vv <- ggvenn(x, fill_color=ggpubr::get_palette(palette="npg", 2),
                 stroke_size = 0.5)
    save_plot(vv, outfile, 5, 5)
}

## Run script
for(feature in c("Gene", "Transcript", "Exon", "Junction")){
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        venn_diagrams(feature, tissue)
    }
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
