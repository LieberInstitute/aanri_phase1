## This script plots enrichment analysis
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg", ".png")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_enrichment <- function(){
    return(data.table::fread("module_enrichment_analysis_DEGs.tsv"))
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(){
    err = 0.0000001
    dt <- load_enrichment() %>% mutate_if(is.character, as.factor) %>%
        mutate(Modules=paste0(Module_ID," (",N_Genes,")"),
               `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) %>%
        mutate(across(Direction, factor,
                      levels=c("All", "Upregulated in AA", "Downregulated in AA")))
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(){
    y0 <- min(memDF()$`log2(OR)`)-0.1
    y1 <- max(memDF()$`log2(OR)`)+0.1
    tile_plot <- memDF() %>%
        ggplot(aes(x = Direction, y = Modules, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        xlab('Ancestry-Related DEGs') + ylab("WGCNA Modules") +
        facet_wrap(.~Tissue, ncol=4, scales="free") + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 15, border=TRUE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold"),
              axis.text.y=element_text(face="bold"),
              strip.text = element_text(face="bold"))
    save_plot(tile_plot, "tileplot_enrichment_modules", 20, 16)
}

## Run script
plot_tile()

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
