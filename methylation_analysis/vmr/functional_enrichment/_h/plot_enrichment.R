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
    err <- 1e-15
    fn  <- paste0(tolower(tissue), "_dmr_rgreat.csv")
    return(data.table::fread(fn) |>
           filter(stringr::str_detect(id, "^GO")) |>
           arrange(p_value) |> head(10) |>
           mutate(`Log10`=-log10(p_adjust+err), Tissue=tissue))
}

generate_dataframe <- function(){
    df_list <- list()
    tissues <- c("Caudate", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list[[jj]] <- get_top_GO(tissues[jj])
    }
    return( bind_rows(df_list) )
}

plot_GO <- function(){
    dt <- generate_dataframe()
    cbPalette <- ggpubr::get_palette(palette = "npg", 4)
    gg1 = ggplot(dt, aes(x=Log10, y=description, color=Tissue,
                         size=fold_enrichment)) +
        geom_point(shape=18, alpha=0.8) +
        labs(y='', x='-Log10 (Adjusted P-value)', size="Enrichment") +
        scale_colour_manual(name="Brain Region", values=cbPalette,
                            labels=c("Caudate","DLPFC","Hippocampus")) +
        scale_size_continuous(range = c(2, 10)) +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

#### MAIN
gg = plot_GO()
save_plot(gg, "ancestry_DMRs.local.GO.stacked", 14, 6)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
