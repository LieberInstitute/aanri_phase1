## Enrichment analysis
library(ggplot2)
library(tidyverse)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_enrichment <- function(feature){
    feat_lt = list("genes"="Gene", "transcripts"="Transcript",
                   "exons"="Exon", "junctions"="Junction")
    fname <- paste0("../../../",feature,
                    "/mainEffect_enrichment/_m/eGene_enrichment_fishers.tsv")
    return(data.table::fread(fname) %>%
           mutate(Feature = feat_lt[[feature]]))
}

gen_data <- function(){
    err = 0.0000001
    features = c("genes", "transcripts", "exons", "junctions")
    dt <- map(features, load_enrichment) %>% bind_rows %>%
        mutate_if(is.character, as.factor) %>%
        mutate(Tissue=fct_relevel(Tissue, rev), `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) %>%
        mutate(across(Feature, factor,
                      levels=c("Gene", "Transcript", "Exon", "Junction")))
    return(dt)
}

plot_tile <- function(dt){
    y0 <- min(dt$`log2(OR)`)-0.1; y1 <- max(dt$`log2(OR)`)+0.1
    tile_plot <- dt %>%
        ggplot(aes(y = Tissue, x = Feature, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        xlab('Ancestry-Related DEGs') + geom_tile(color = "grey") +
        ylab('') + ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(-y1, 0, y1)),
                             limits=c(-y1,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold"),
              axis.text=element_text(face="bold"),
              strip.text = element_text(face="bold"))
    save_plot(tile_plot, "tileplot_enrichment_eGenes", 8, 6)
}

## Run script
dt <- gen_data(); plot_tile(dt)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
