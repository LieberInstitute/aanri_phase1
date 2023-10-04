#### This script plots the distribution of AFD

library(ggpubr)

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_data <- function(tissue, feature){
    new_tissue <- gsub(" ", "_", tolower(tissue))
    fn <- paste(new_tissue,feature,"degs_AFD.tsv", sep=".")
    return(data.table::fread(fn) |>
           dplyr::select(gene_id, AFD, DEG) |>
           dplyr::mutate(region=tissue))
}

merge_regions <- function(feature){
    regions = c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
    purrr::map2(regions, feature, \(x, y) load_data(x, y))
    return(purrr::map2(regions, feature, \(x, y) load_data(x, y)) |>
           dplyr::bind_rows())
}

plot_dist <- function(feature){
    outfile <- paste0("density_plot.allele_freq.", feature)
    df <- merge_regions(feature) |>
        dplyr::mutate(DE=ifelse(DEG == 1, "DEG", "Other"))
    dens <- ggdensity(df, x="AFD", fill="DE", color="DE", facet.by="region",
                      add="mean", rug=FALSE, ylab="", palette="npg",
                      panel.labs.font=list(face="bold"),
                      xlab="Allele Frequency Difference",
                      ggtheme=theme_pubr(base_size=15), ncol=4)
    save_plot(dens, outfile, 12, 4)
}


#### Main
plot_dist("eGene")
plot_dist("eSNP")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
