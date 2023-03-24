## This script generates plots for GSEA analysis.
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

get_global_dmrs <- function(tissue){
    fn <- paste0("../../_m/aaOnly/",
                 tolower(tissue), "_dmr.csv")
    return(data.table::fread(fn) |>
           mutate(id=paste(seqnames, start, end, sep="_")) |>
           select(id, beta, fdr))
}

get_local_dmrs <- function(tissue){
    fn <- paste0("../../_m/aaOnly/",
                 tolower(tissue), "_dmr_local.csv")
    return(data.table::fread(fn) |>
           mutate(id=paste(seqnames, start, end, sep="_")) |>
           select(id, beta, fdr))
}

merge_data <- function(tissue){
    return(inner_join(get_global_dmrs(tissue),
                      get_local_dmrs(tissue), by="id",
                      suffix=c("_global", "_local")) |>
           mutate(region = tissue))
}

generate_dataframe <- function(SIG){
    df_list <- list()
    tissues <- c("Caudate", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list[[jj]] <- merge_data(tissues[jj])
    }
    if(SIG){
        return( bind_rows(df_list) |>
               filter(fdr_global < 0.05 | fdr_local < 0.05) )
    } else {
        return( bind_rows(df_list) )
    }
}

plot_correlation <- function(SIG=FALSE){
    dt  <- generate_dataframe(SIG)
    sca <- ggscatter(dt, x="beta_global", y="beta_local",
                     color="black", alpha=0.45,
                     panel.labs.font=list(face="bold", size=16),
                     facet.by="region", add="reg.line",
                     xlab="Global Ancestry (Effect Size)",
                     ylab="Local Ancestry (Effect Size)",
                     add.params=list(color="blue", fill="lightgray"),
                     conf.int=TRUE, cor.coef=TRUE, ncol=3,
                     ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=16)
    return(sca)
}

#### MAIN
                                        # All DMRs
gg <- plot_correlation()
save_plot(gg, "ancestry_DMRs.comparison.scatter.all", 14, 5)
                                        # Significant DMRs
gg <- plot_correlation(TRUE)
save_plot(gg, "ancestry_DMRs.comparison.scatter.signif", 14, 5)

                                        # Save results
generate_dataframe() |>
    data.table::fwrite("dmr_global_local_combined_3regions.csv")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
