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

get_top_GO <- function(tissue, mhc_region){
    new_tissue <- gsub(" ", "_", tolower(tissue))
    filenames <- list.files(path=paste0("../_m/", mhc_region, "/"),
                            full.names=TRUE,
                            pattern=paste0(new_tissue, "_enrich*"))
    filenames <- filenames[!grepl("CC", filenames)]
    df_list <- list()
    for(i in seq_along(filenames)){
        df_list[[i]] <- data.table::fread(filenames[i])
    }
    return(bind_rows(df_list) |> arrange(pvalue) |> head(10) |>
           mutate(`Log10`=-log10(qvalue), Tissue=tissue))
}

generate_dataframe <- function(mhc_region){
    df_list <- list()
    tissues <- c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list[[jj]] <- get_top_GO(tissues[jj], mhc_region)
    }
    return( bind_rows(df_list) )
}

plot_GO <- function(mhc_region){
    dt <- generate_dataframe(mhc_region)
    cbPalette <- ggpubr::get_palette(palette = "npg", 4)
    gg1 = ggplot(dt, aes(x=Log10, y=Description, color=Tissue)) +
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
tissues <- c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")
for(mhc_region in c("HLA", "MHC", "xMHC")){
                                        # GO plotting
    gg = plot_GO(tolower(mhc_region))
    fn = paste0("ancestry_mash_GOenrich_top10_stacked.no", mhc_region)
    save_plot(gg, fn, 8, 6)
                                        # Similarity
    outfile <- paste0("GO_similarity_summary.no",mhc_region,".tsv")
    purrr::map(tissues, similarity_go) |> bind_rows() |>
        data.table::fwrite(outfile, sep="\t")

}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
