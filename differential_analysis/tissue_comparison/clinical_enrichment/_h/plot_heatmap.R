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

plot_tile <- function(label){
    tile_plot <- memENRICH() %>% mutate_if(is.character, as.factor) %>%
        filter(str_detect(Comparison, label)) %>%
        mutate(Tissue=fct_relevel(Tissue, rev), `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) %>%
        ggplot(aes(x = Comparison, y = Tissue, fill = `OR Percentile`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        ylab('Ancestry-Related DEGs') + xlab('') +
        geom_tile(color = "grey") + ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradient2(midpoint = 0.5, low = "blue", mid = "white",
                             high = "red", space = "Lab", limits=c(0,1)) +
        ggpubr::theme_pubr(base_size = 20, border=TRUE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold"),
              axis.text.y=element_text(face="bold"))
    save_plot(tile_plot, paste0("tileplot_enrichment_",tolower(label)), 10, 7)
}

for(label in c("DEG", "TWAS")){
    plot_tile(label)
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
