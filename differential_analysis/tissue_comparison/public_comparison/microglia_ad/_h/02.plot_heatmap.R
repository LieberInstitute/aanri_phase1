#### Plotting heatmap

library(ggplot2)
library(dplyr)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_enrichment <- function(){
    return(data.table::fread("celltype_enrichment.DEGs_states.tsv"))
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(){
    err = 0.0000001
    dt <- memENRICH() |> 
        mutate_if(is.character, as.factor) |>
        mutate(`-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE))
    mstates <- paste0("MG", seq(0,12))[!grepl("MG9", paste0("MG", seq(0,12)))]
    dt$State <- factor(dt$State, levels=mstates)
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(label, w, h){
    y0 <- min(memDF()$`log2(OR)`)-0.1
    y1 <- max(memDF()$`log2(OR)`)+0.1
    tile_plot <- memDF() |>
        mutate(State=forcats::fct_relevel(State, rev)) |>
        ggplot(aes(y = State, x = Tissue, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        xlab('Ancestry-Associated DEGs') + ylab("Microglial States") +
        facet_wrap(Response~Direction) + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=TRUE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              axis.title=element_text(face="bold"),
              strip.text=element_text(face="bold"))
    save_plot(tile_plot, paste0("tileplot_enrichment.",tolower(label)), w, h)
}

## Run script
plot_tile("adDEG_microglia", 10, 12)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
