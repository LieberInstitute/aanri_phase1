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

load_data <- function(feature){
    if(tolower(feature) == "gene"){
        return(data.table::fread("brainseq_degs_constrain_score.tsv"))
    } else {
        return(data.table::fread("brainseq_degs_constrain_score_tx.tsv"))
    }
}

plot_errplot <- function(feature){
    xlabel = paste0("LOEUF decile bin (",feature,")")
    ylabel = ifelse(tolower(feature) == "gene",
                    "Ancestry-Associated DEGs\n(Mean lfsr)",
                    "Ancestry-Associated DE\n(Mean lfsr)")
    errplt <- load_data(feature) %>%
        #mutate(LOEUF_decile = (oe_lof_upper_bin + 0.5) * 10) %>%
        ggerrorplot(x="oe_lof_upper_bin", y="lfsr", color="Tissue",
                    palette="npg", ylab=ylabel, xlab=xlabel,
                    desc_stat="mean_ci", error.plot="pointrange",
                    position=position_dodge(0.5),
                    ggtheme=theme_pubr(base_size=15)) +
        font("xy.title", face="bold", size=20)
    save_plot(errplt, paste0("constrain_correlation_errorplot_",
                             tolower(feature)), 7, 7)
}

plot_corr <- function(feature){
    ylabel = paste(feature, "Constrain (LOEUF)")
    xlabel = ifelse(tolower(feature) == "gene", "DEGs (lfsr)", "DE (lfsr)")
    corr <- load_data(feature) %>%
        ggscatter(x="lfsr", y="oe_lof_upper", facet.by="Tissue",
                  color="black", alpha=0.2, xlab=xlabel, ylab=ylabel,
                  add="reg.line", panel.labs.font=list(face="bold"),
                  add.params=list(color="blue", fill="lightgray"),
                  ylim=c(-2,4), conf.int=TRUE, cor.coef=TRUE,
                  cor.coeff.args=list(method="pearson", label.x=0.10,
                                      label.y=3, label.sep="\n"),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=4) +
        font("xy.title", face="bold", size=20)
    save_plot(corr, paste0("constrain_correlation_lfsr_",tolower(feature)), 16, 4)
}

plot_density <- function(feature, variable){
    xlabel = ifelse(tolower(feature) == "gene", "DEGs (lfsr)", "DE (lfsr)")
    fn     = paste("constrain_density_lfsr",tolower(variable),
               tolower(feature),sep="_")
    den <- load_data(feature) %>% mutate_at("oe_lof_upper_bin",as.character) %>%
        mutate("Constrained"=ifelse(oe_lof_upper > 0.75, "No", "Yes")) %>%
        ggdensity(x="lfsr", color=variable, fill=variable, facet.by="Tissue",
                  xlab=xlabel, palette="npg", rug=TRUE, add="mean",
                  panel.labs.font=list(face="bold"), ncol=4, 
                  ggtheme=theme_pubr(base_size=15, border=TRUE)) +    
        font("xy.title", face="bold", size=20)
    save_plot(den, fn, 16, 4)
}

#### MAIN
plot_errplot("Gene")
plot_errplot("Transcript")
plot_corr("Gene")
plot_corr("Transcript")
plot_density("Gene", "Constrained")
plot_density("Gene", "oe_lof_upper_bin")
plot_density("Transcript", "Constrained")
plot_density("Transcript", "oe_lof_upper_bin")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
