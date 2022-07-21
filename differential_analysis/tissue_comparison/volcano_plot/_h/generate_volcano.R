## Generating fancy volcano plot
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

get_mash_lfsr <- function(tissue, feature){
    fn <- paste0("../../_m/",feature,"/lfsr_feature_4tissues.txt.gz")
    return(data.table::fread(fn) %>% select(all_of(c("Effect", tissue))) %>%
           dplyr::rename("lfsr"=tissue))
}

get_de <- function(tissue, feature){
    fn <- paste0("../../../",tissue,"/_m/",feature,
                 "/diffExpr_EAvsAA_full.txt")
    return(data.table::fread(fn) %>% janitor::clean_names() %>%
           mutate(fold_change = 2^log_fc) %>%
           select(v1,fold_change,adj_p_val))
}

get_data <- function(tissue, feature){
    new_tissue <- ifelse(tissue=="Dentate.Gyrus","dentateGyrus",tolower(tissue))
    return(get_de(new_tissue, feature) %>%
           inner_join(get_mash_lfsr(tissue, feature),by=c("v1"="Effect")) %>%
           dplyr::rename("feature_id"="v1") %>%
           mutate(gene_type = case_when(fold_change >= 2 & lfsr <= 0.05 ~ "up",
                                        fold_change <= 0.5 & lfsr <=0.05 ~ "down",
                                        TRUE ~ "ns")))
}

plot_volcano <- function(tissue, feature){
    cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
    vol_plot <- get_data(tissue, feature) %>%
        ggplot(aes(x=log2(fold_change), y=-log10(lfsr))) +
        geom_point(aes(color=gene_type), alpha=0.3,shape=16,size=1) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed") +
        geom_vline(xintercept = c(log2(0.5), log2(2)), linetype="dashed") +
        scale_color_manual(values=cols) +
        scale_x_continuous(breaks=c(seq(-10,10,2)), limits=c(-10, 10)) +
        labs(title=toupper(feature), x="log2(fold change)",
             y="-log10(lfsr)", color="Expression change") +
        theme_bw(base_size=15) +
        theme(plot.title = element_text(hjust=0.5),
              panel.border = element_rect(color="black",fill=NA,size=0.5),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    return(vol_plot)
}

figure_annotation <- function(tissue){
    gg <- plot_volcano(tissue, "genes")
    tt <- plot_volcano(tissue, "transcripts")
    ee <- plot_volcano(tissue, "exons")
    jj <- plot_volcano(tissue, "junctions")
    figure <- ggarrange(gg, tt, ee, jj, ncol=4, nrow=1,common.legend=TRUE,
                        font.label=list(face="bold"))
    new_fig <- annotate_figure(figure,
                               left=text_grob(tissue,rot=90,size=20,face="bold"))
    return(new_fig)
}

#### MAIN
for(tissue in c("Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus")){
    new_fig <- figure_annotation(tissue)
    ggsave(new_fig, filename=paste(tolower(tissue), "volcano_plot.pdf", sep="_"),
           width=16, height=4)
    ggsave(new_fig, filename=paste(tolower(tissue), "volcano_plot.png", sep="_"),
           width=16, height=4)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
