#### Generate enrichment plot
library(ggplot2)
library(dplyr)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_enrichment <- function(label){
    fn <- paste0("xmhc_enrichment_analysis.",label,".txt")
    return(data.table::fread(fn))
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(label){
    err = 0.0000001
    dt <- memENRICH(label) |> mutate_if(is.character, as.factor) |>
        mutate(OR=ifelse(is.na(OR), 1, OR),
               OR=ifelse(is.infinite(OR), 30, OR)) |>
        mutate(`-log10(FDR)`= -log10(FDR),
               `-log10(p-value)`= -log10(PValue),
               `OR Percentile`= OR / (1+OR),
               p.fdr.sig=FDR < 0.05, p.sig=PValue < 0.01,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) |>
        mutate(across(Direction, \(x) factor(x, levels=c("All", "Down", "Up"))))
    levels(dt$Direction) <- c("All", "AA Bias", "EA Bias")
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(label, w, h){
    y0 <- min(memDF(tolower(label))$`log2(OR)`)-0.1
    y1 <- max(memDF(tolower(label))$`log2(OR)`)+0.1
    tile_plot <- memDF(tolower(label)) |>
        ggplot(aes(y = .data[[label]], x = Region, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        xlab('Ancestry-Associated DEGs') + facet_grid(.~Direction) +
        geom_tile(color = "grey") + ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold", size=22),
              strip.text=element_text(face="bold", size=22))
    save_plot(tile_plot, paste0("tileplot_enrichment_xmhc.",tolower(label)), w, h)
}

plot_tile_pvalue <- function(label, w, h){
    y0 <- min(memDF(tolower(label))$`log2(OR)`)-0.1
    y1 <- max(memDF(tolower(label))$`log2(OR)`)+0.1
    tile_plot <- memDF(tolower(label)) |>
        ggplot(aes(y = .data[[label]], x = Region, fill = `log2(OR)`,
                   label = ifelse(p.sig,
                                  format(round(`-log10(p-value)`,1), nsmall=1), ""))) +
        xlab('Ancestry-Associated DEGs') + facet_grid(.~Direction) +
        geom_tile(color = "grey") + ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold", size=22),
              strip.text=element_text(face="bold", size=22))
    save_plot(tile_plot, paste0("tileplot_enrichment_xmhc.",tolower(label)), w, h)
}

## Run script
plot_tile("Subregion", 15, 6)
plot_tile_pvalue("Gene_Cluster", 15, 8)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
