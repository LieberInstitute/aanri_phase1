## Heatmap of enrichment
suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyverse)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_enrichment <- function(GENE){
    if(GENE){
        return(data.table::fread("constrain_enrichment.tsv"))
    } else {
        dt <- data.table::fread("constrain_enrichment_tx.tsv") %>%
            mutate_if(is.character, as.factor) %>%
            mutate(Canonical=factor(Canonical, levels=c("All", "True", "False")))
        levels(dt$Canonical) <- c("All", "Canonical", "Non-Canonical")
        return(dt)
    }
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(GENE){
    err = 0.0000001
    dt <- memENRICH(GENE) %>% mutate_if(is.character, as.factor) %>%
        mutate_at("Upper_Bin", as.factor) %>%
        mutate(Upper_Bin=fct_relevel(Upper_Bin, rev), `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= Odds_Ratio / (1+Odds_Ratio), p.fdr.sig=FDR < 0.05,
               `Log2(OR)` = log2(Odds_Ratio+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE))
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(w, h, GENE=TRUE){
    y0 <- min(memDF(GENE)$`Log2(OR)`)-0.1; y1 <- max(memDF(GENE)$`Log2(OR)`)+0.1
    tile_plot0 <- memDF(GENE) %>%
        ggplot(aes(y=Upper_Bin, x=Tissue, fill=`Log2(OR)`,
                   label = ifelse(p.fdr.sig,format(round(`-log10(FDR)`,1),nsmall=1),""))) +
        xlab('Ancestry-Related DEGs') + ylab("Constrain (LOEUF)")
    if(GENE){
        tile_plot0 <- tile_plot0
        label <- "genes"
    } else {
        tile_plot0 <- tile_plot0 + facet_grid(.~Canonical)
        label <- "tx"
    }
    tile_plot <- tile_plot0 + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, abs(y0))),
                             limits=c(y0,abs(y0))) +
        ggpubr::theme_pubr(base_size = 20, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold", size=22),
              strip.text=element_text(face="bold", size=22))
    save_plot(tile_plot, paste0("tileplot_enrichment_",label), w, h)
}

## Run script
plot_tile(6, 8, TRUE)
plot_tile(10, 8, FALSE)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
