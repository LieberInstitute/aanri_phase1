## Generate plots for LDSC results
library(ggplot2)
library(tidyverse)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_filtered_ldsc <- function(){
    return(data.table::fread("filtered_ldsc_results.tsv"))
}
memLDSC <- memoise::memoise(load_filtered_ldsc)

gen_data <- function(feature){
    dt <- memLDSC() %>% filter(Feature == feature) %>%
        mutate_if(is.character, as.factor) %>%
        mutate(trait=fct_relevel(trait, rev),
               `-log10(FDR)`= -log10(Enrichment_FDR),
               p.fdr.sig=Enrichment_FDR < 0.05,
               p.fdr.cat=cut(Enrichment_FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) %>%
        mutate(across(Direction, factor, levels=c("", "down", "up")))
    levels(dt$Direction) <- c("All", "Increased in AA", "Decreased in AA")
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(feature, w, h){
    if(feature == "Gene"){
        xlabel = 'Ancestry-Related DEGs'
    } else {
        xlabel = 'Ancestry-Related DETs'
    }
    ##y0 <- min(memDF(feature)$Enrichment)-0.1
    y0 <- -0.25
    y1 <- max(memDF(feature)$Enrichment)+0.1
    tile_plot <- memDF(feature) %>%
        ggplot(aes(y = trait, x = Tissue, fill = Enrichment,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1),
                                         nsmall=1), ""))) +
        xlab(xlabel) + ylab("") + facet_grid(.~Direction) +
        geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0, y1)) +
        ggpubr::theme_pubr(base_size = 15, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold", size=18),
              strip.text=element_text(face="bold", size=16))
    save_plot(tile_plot, paste0("tileplot_ldsc_enrichment_",tolower(feature)), w, h)
}

## Run script
plot_tile("Gene", 12, 8)
plot_tile("Transcript", 12, 8)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
