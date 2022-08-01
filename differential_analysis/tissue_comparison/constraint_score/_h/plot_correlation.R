## This script plots the correlation between lfsr and constrain score
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    return(data.table::fread("brainseq_degs_constrain_score.tsv"))
}

plot_corr <- function(){
    corr <- load_data() %>%
        ggscatter(x="lfsr", y="oe_lof_upper", facet.by="Tissue",
                  color="black", alpha=0.2, xlab="DEGs (lfsr)",
                  ylab="Gene Constrain (LOEUF)", add="reg.line",
                  panel.labs.font=list(face="bold"), ylim=c(-2,4),
                  add.params=list(color="blue", fill="lightgray"),
                  conf.int=TRUE, cor.coef=TRUE,
                  cor.coeff.args=list(method="pearson", label.x=0.10,
                                      label.y=3, label.sep="\n"),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=4) +
        font("xy.title", face="bold", size=20)
    save_plot(corr, "constrain_correlation_lfsr", 16, 4)
}

#### MAIN
plot_corr()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
