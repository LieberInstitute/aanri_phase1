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


get_top_GO <- function(tissue, SHARED=FALSE, ANCESTRY=TRUE){
    new_tissue <- gsub(" ", "_", tolower(tissue))
    pat <- ifelse(SHARED, paste0(new_tissue, "_.*shared.tsv$"),
                  paste0(new_tissue, "_.*unique.tsv$"))
    filenames <- ifelse(ANCESTRY,
                        list.files(path="../_m/ancestry", full.names=TRUE, pattern=pat),
                        list.files(path="../_m/environmental", full.names=TRUE, pattern=pat))
    df_list <- list()
    for(i in seq_along(filenames)){
        df_list[[i]] <- data.table::fread(filenames[i])
    }
    return(bind_rows(df_list) |> arrange(pvalue) |> head(10) |>
           mutate(`Log10`=-log10(qvalue), Tissue=tissue))
}

generate_dataframe <- function(SHARED){
    df_list <- list()
    tissues <- c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
    for(ii in seq_along(tissues)){
        dt_list <- list(); ancestry <- c(TRUE, FALSE)
        for(jj in seq_along(ancestry)){
            dt_list[[jj]] <- get_top_GO(tissues[ii], SHARED, ancestry[jj]) |>
                mutate(deg_source=ifelse(ancestry[jj], "Ancestry", "Environmental"))
        }
        df_list[[ii]] <- bind_rows(dt_list)
    }
    return( bind_rows(df_list) )
}

plot_GO <- function(SHARED=FALSE){
    dt <- generate_dataframe(SHARED)
    cbPalette <- ggpubr::get_palette(palette = "npg", 4)
    gg1 = ggplot(dt, aes(x=Log10, y=Description, color=Tissue)) +
        facet_wrap("~deg_source") + 
        geom_point(shape=18, alpha=0.8, size=5) + #xlim(-2, 2) +
        labs(y='', x='-Log10 (q-value)') +
        scale_colour_manual(name="Brain Region", values=cbPalette,
                            labels=c("Caudate", "Dentate Gyrus",
                                     "DLPFC", "Hippocampus")) +
        geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

#### MAIN
                                        # Unique
gg = plot_GO(FALSE)
save_plot(gg, "ancestry_mash_GOenrich_top10_stacked.unique", 12, 9)
                                        # Shared
gg = plot_GO(TRUE)
save_plot(gg, "ancestry_mash_GOenrich_top10_stacked.shared", 16, 9)


#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
