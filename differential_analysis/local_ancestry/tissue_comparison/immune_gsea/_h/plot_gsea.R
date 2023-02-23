## This script generates plots for GSEA analysis.
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}


get_top_GO <- function(tissue){
    new_tissue <- gsub(" ", "_", tolower(tissue))
    filenames <- list.files(path="../_m/genes/", full.names=TRUE,
                            pattern=paste0(new_tissue, "_gsea*"))
    df_list <- list()
    for(i in seq_along(filenames)){
        df_list[[i]] <- data.table::fread(filenames[i])
    }
    return(bind_rows(df_list) %>% arrange(pvalue) %>%
           filter(stringr::str_detect(Description, "immune | T cell | B cell")) %>%
           mutate(`-Log10(qvalues)`=-log10(qvalues), Tissue=tissue))
}

generate_dataframe <- function(){
    df_list <- list()
    tissues <- c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list[[jj]] <- get_top_GO(tissues[jj])
    }
    go_terms <- bind_rows(df_list) %>% filter(qvalues < 0.25) %>% pull(ID) %>% unique
    return( bind_rows(df_list) %>% filter(ID %in% go_terms) %>%
            mutate(p.fdr.sig=qvalues < 0.25,
                   p.fdr.cat=cut(qvalues, breaks=c(1,0.05,0.01,0.005,0),
                                 labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                                 include.lowest=TRUE)) )
}

plot_tile <- function(dt){
    y0 <- min(dt$NES)-0.1; y1 <- max(dt$NES)+0.1
    tile_plot <- dt %>%
        ggplot(aes(x=Tissue, y=Description, fill=NES,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-Log10(qvalues)`,1), nsmall=1), ""))) +
        xlab('Ancestry-Related DEGs') + ylab('') + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 15, border=TRUE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold", size=18),
              strip.text=element_text(face="bold", size=16))
    return(tile_plot)
}

#### MAIN
                                        # Immune related data
dt <- generate_dataframe()
dt %>% data.table::fwrite("immune_related_GSEA_results.tsv", sep='\t')
                                        # Tile plot
tile_plot <- plot_tile(dt)
save_plot(tile_plot, "tileplot_immune_GSEA_enrichment", 16, 6)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
