## This script generates plots of the metrics

library(ggpubr)

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

#### Main

df <- data.table::fread("partial_r2_de_summary.tsv") |>
    dplyr::filter(feature_type == "Gene")

hist <- gghistogram(df, x="delta_partial", fill="lightgray", add="mean",
                    rug=TRUE, add_density=TRUE, facet.by="region",
                    panel.labs.font=list(face="bold"), 
                    palette="npg", xlab="Partial R2", ylab="Frequencey",
                    ggtheme=theme_pubr(base_size=15, border=TRUE))
save_plot(hist, "distribution_partial_r2", 10, 4)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
