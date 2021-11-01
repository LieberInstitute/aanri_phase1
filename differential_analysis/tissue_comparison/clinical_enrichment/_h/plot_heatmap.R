library(ggplot2)
library(tidyverse)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_enrichment <- function(){
    return(data.table::fread("clincial_phenotypes_enrichment_analysis.tsv"))
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(label){
    err = 0.0000001
    dt <- memENRICH() %>% mutate_if(is.character, as.factor) %>%
        filter(str_detect(Comparison, label)) %>%
        mutate(Tissue=fct_relevel(Tissue, rev), `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) %>%
        mutate(across(Direction, factor, levels=c("All", "AA Bias", "EA Bias")))
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(label){
    y0 <- min(memDF(label)$`log2(OR)`)-0.1
    y1 <- max(memDF(label)$`log2(OR)`)+0.1
    tile_plot <- memDF(label) %>%
        ggplot(aes(x = Comparison, y = Tissue, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        ylab('Ancestry-Related DEGs') + xlab('') + facet_grid(.~Direction) +
        geom_tile(color = "grey") + ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=TRUE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold"),
              axis.text.y=element_text(face="bold"))
    save_plot(tile_plot, paste0("tileplot_enrichment_",tolower(label)), 14, 7)
}

for(label in c("DEG", "TWAS")){
    plot_tile(label)
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
