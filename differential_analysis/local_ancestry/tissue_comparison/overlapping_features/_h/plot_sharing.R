## Generate plots comparing local and global ancestry results
##  1. heatmap of enrichment
##  2. venn diagram showing overlap

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

get_local_mash <- function(feature, tissue){
    fn <- paste0("../../summary_table/_m/",
                 "BrainSeq_ancestry_4features_4regions.txt.gz")
    return(data.table::fread(fn) |>
           filter(feature_type == feature, region == tissue) |>
           select(feature_id, posterior_mean, lfsr))
}
memLOCAL <- memoise::memoise(get_local_mash)
                 
get_global_mash <- function(feature, tissue){
    fn <- here::here("differential_analysis/tissue_comparison",
                     "summary_table/_m",
                     "BrainSeq_ancestry_4features_4regions.txt.gz")
    return(data.table::fread(fn) |>
           filter(Type == feature, Tissue == tissue) |>
           select(Effect, posterior_mean, lfsr) |>
           rename("feature_id"="Effect"))
}
memGLOBAL <- memoise::memoise(get_global_mash)

load_enrichment <- function(){
    return(data.table::fread("enrichment_analysis.txt"))
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(){
    err = 0.0000001
    dt <- memENRICH() |> mutate_if(is.character, as.factor) |>
        mutate(Feature=forcats::fct_relevel(Feature, rev),
               `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) |>
        mutate(across(Direction, \(x) factor(x, levels=c("All", "Down", "Up"))))
    levels(dt$Direction) <- c("All", "AA Bias", "EA Bias")
    return(dt)
}
memDF <- memoise::memoise(gen_data)

venn_diagrams <- function(feature, tissue){
    new_feature <- tolower(feature)
    new_tissue  <- gsub(" ", "_", tolower(tissue))
    outfile     <- paste("venn_diagram", new_feature, new_tissue, sep="_")
    x  <- list(
        LOCAL  = memLOCAL(feature, tissue)$feature_id,
        GLOBAL = memGLOBAL(feature, tissue)$feature_id
    )
    vv <- ggvenn(x, fill_color=ggpubr::get_palette(palette="npg", 2),
                 stroke_size = 0.5)
    save_plot(vv, outfile, 5, 5)
}

plot_tile <- function(w, h){
    y0 <- min(memDF()$`log2(OR)`)-0.1
    y1 <- max(memDF()$`log2(OR)`)+0.1
    tile_plot <- memDF() |>
        ggplot(aes(x = Feature, y = Region, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        ylab('Ancestry-Related DEGs') + xlab("") + facet_grid(.~Direction) +
        geom_tile(color = "grey") + ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(-y1, 0, y1)),
                             limits=c(-y1,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold"),
              strip.text=element_text(face="bold"))
    save_plot(tile_plot, "tileplot_enrichment_ancestry", w, h)
}

## Run script
for(feature in c("Gene", "Transcript", "Exon", "Junction")){
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        venn_diagrams(feature, tissue)
    }
}

plot_tile(10, 6)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
