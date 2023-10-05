## Generate plots for LDSC results
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_ldsc <- function(){
    traits <- c("Alzheimer Disease3", "Parkinson Disease")
    return(data.table::fread("../_h/ldsc_results.txt") |>
           mutate(FDR = p.adjust(Enrichment_p, method="fdr")) |>
           filter(trait %in% traits) |> as.data.frame() |>
           mutate(cell_type=ifelse(grepl("^MG", group),
                                   "Microglia",
                            ifelse(grepl("^AST", group),
                                   "Astrocyte", "Oligodendrocyte"))))
}
memLDSC <- memoise::memoise(load_ldsc)

gen_data <- function(){
    dt <- memLDSC() |> mutate_if(is.character, as.factor) |>
        mutate(trait=forcats::fct_relevel(trait, rev),
               `-log10(FDR)`= -log10(FDR),
               p.fdr.sig=FDR < 0.01,
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE))
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(w, h){
    xlabel  <- 'Glia Subpopulation'
    outfile <- "tileplot.ldsc_enrichment.cell_type"
    y1 <- max(memDF()$Enrichment)+0.1
    tile_plot <- memDF() |>
        ggplot(aes(y = trait, x = group, fill = Enrichment,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1),
                                         nsmall=1), ""))) +
        xlab(xlabel) + ylab("") + geom_tile(color = "grey") +
        facet_wrap(.~cell_type, scales="free", ncol=1) +
        ggfittext::geom_fit_text(contrast=TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(0, 1, y1)),
                             limits=c(0, y1)) +
        ggpubr::theme_pubr(base_size = 15, border=TRUE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold", size=18),
              strip.text=element_text(face="bold", size=16))
    save_plot(tile_plot, outfile, w, h)
}

## Main
plot_tile(7, 8)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
