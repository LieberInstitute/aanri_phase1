## Plotting enrichment
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_enrichment <- function(){
    return(data.table::fread("nedelec_immune_enrichment_analysis_DEGs.tsv"))
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(){
    err = 0.0000001
    dt <- memENRICH() %>% mutate_if(is.character, as.factor) %>%
        mutate(`-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) %>%
        mutate(across(Direction, factor, levels=c("All", "Down", "Up")))
    levels(dt$Direction) <- c("All", "Upregulated in AA", "Downregulated in AA")
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(){
    y1 <- max(memDF()$`log2(OR)`)+0.1
    tile_plot <- memDF() %>%
        mutate(Conditions = factor(Conditions,
                                   levels=c("Non-infected", "Listeria", "Salmonella"))) %>%
        ggplot(aes(y = Conditions, x = Tissue, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        xlab('Ancestry-Related DEGs') + ylab("Nédélec et al. pop-DE") +
        facet_grid(.~Direction) + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(-y1/2, 0, y1)),
                             limits=c(-y1/2,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold", size=24),
              strip.text = element_text(face="bold"))
    return(tile_plot)
}

## MAIN
tile_plot <- plot_tile()
save_plot(tile_plot, "tileplot_enrichment_nedelec", 12, 5)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
