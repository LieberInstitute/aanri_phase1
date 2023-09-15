## This script generates plots for GSEA analysis.
suppressPackageStartupMessages({
    library(dplyr)
    library(rrvgo)
    library(ggplot2)
})

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

similarity_go <- function(tissue){
    new_tissue <- gsub(" ", "_", tolower(tissue))
    df_list <- list()
    for(ONT in c("BP", "MF", "CC")){
        fn <- paste0("genes/", new_tissue,"_enrich_",ONT,".tsv")
        df <- data.table::fread(fn)
        simMat <- calculateSimMatrix(df$ID, orgdb="org.Hs.eg.db",
                                     ont=ONT, method="Rel")
        scores <- setNames(-log10(df$qvalue), df$ID)
        reTerm <- reduceSimMatrix(simMat, scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
        reTerm["Ont"]  <- ONT
        df_list[[ONT]] <- reTerm
    }
    return(bind_rows(df_list) |> mutate(Tissue=tissue))
}

get_top_GO <- function(tissue){
    new_tissue <- gsub(" ", "_", tolower(tissue))
    filenames <- list.files(path="../_m/genes/", full.names=TRUE,
                            pattern=paste0(new_tissue, "_enrich*"))
    filenames <- filenames[!grepl("CC", filenames)]
    df_list <- list()
    for(i in seq_along(filenames)){
        df_list[[i]] <- data.table::fread(filenames[i])
    }
    return(bind_rows(df_list) |> arrange(pvalue) |> head(10) |>
           mutate(`Log10`=-log10(qvalue), Tissue=tissue))
}

generate_dataframe <- function(){
    df_list <- list()
    tissues <- c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list[[jj]] <- get_top_GO(tissues[jj])
    }
    return( bind_rows(df_list) )
}

plot_GO <- function(){
    dt <- generate_dataframe()
    cbPalette <- ggpubr::get_palette(palette = "npg", 4)
    gg1 = ggplot(dt, aes(x=Log10, y=Description, color=Tissue)) +
        geom_point(shape=18, alpha=0.8, size=5) + #xlim(-2, 2) +
        labs(y='', x='-Log10 (q-value)') +
        scale_colour_manual(name="Brain Region", values=cbPalette,
                            labels=c("Caudate", "Dentate Gyrus",
                                     "DLPFC", "Hippocampus")) +
        ## scale_size_continuous(range = c(2, 10)) +
        geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

#### MAIN
gg = plot_GO()
save_plot(gg, "ancestry_mash_GOenrich_top10_stacked", 9, 6)

                                        # Similarity
tissues <- c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
purrr::map(tissues, similarity_go) |> bind_rows() |>
    data.table::fwrite("GO_similarity_summary.tsv", sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
