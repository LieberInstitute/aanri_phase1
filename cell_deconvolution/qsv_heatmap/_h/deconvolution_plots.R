## Author: Louise Huuki
## Edited by KJ Benjamin

library(SummarizedExperiment)
library(DeconvoBuddies)
library(RColorBrewer)
library(sessioninfo)
library(patchwork)
library(tidyverse)
library(viridis)
library(broom)

## Functions
get_qsv_file <- function(tissue){
    qsv_lt = list("Caudate"="../../../input/phenotypes/_m/qSV_caudate.csv",
                  "DLPFC"="../../../input/phenotypes/_m/qSV_dlpfc.csv",
                  "Hippocampus"="../../../input/phenotypes/_m/qSV_hippo.csv",
                  "Dentate Gyrus"="../../../input/phenotypes/_m/qSV_dg.csv")
    return(qsv_lt[[tissue]])
}

save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_pheno <- function(){
    fname = "../../../input/phenotypes/merged/_m/merged_phenotypes.csv"
    df = data.table::fread(fname) %>% select(-V1) %>%
        filter(Dx %in% c("Control"), Age > 17, Race %in% c("AA", "CAUC")) %>%
        mutate(Race = gsub("CAUC", "EA", Race))
    return(df)
}
memPHENO <- memoise::memoise(get_pheno)

get_cell_prop <- function(){
    ## Load Bisque Estimated Props
    load("../../_m/est_prop_Bisque.Rdata")
    cc = est_prop_bisque$caudate$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="Caudate")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="DLPFC")
    hh = est_prop_bisque$hippo$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="Hippocampus")
    gg = est_prop_bisque$dg$Est.prop.long %>%
        separate(sample, c("sample", "batch")) %>%
        inner_join(memPHENO(), by=c("sample"="RNum")) %>%
        mutate(Tissue="Dentate Gyrus")
    df = bind_rows(cc, dd, hh, gg)
    return(df)
}
memPROP <- memoise::memoise(get_cell_prop)

#### Load colors and plotting functions ####
cell_colors <- create_cell_colors(pallet = "classic")
cell_colors <- cell_colors[c("Astro", "Endo", "Micro", "Mural", "Oligo",
                             "OPC", "Tcell", "Excit", "Inhib")]
#### Comparison plots ####
for(tissue in c("Caudate", "DLPFC", "Hippocampus", "Dentate Gyrus")){
    dirname = gsub(" ", "_", tolower(tissue))
    dir.create(dirname)
    ## Boxplot
    group_boxPlot <- memPROP() %>% filter(Tissue == tissue) %>%
        ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
        geom_boxplot() + labs(x = "Cell Type", y = "Proportion") +
        ggpubr::theme_pubr(base_size = 15) +
        scale_fill_manual(values = cell_colors, guide = "none")
    ## Save image
    save_img(group_boxPlot, paste0(dirname,"/cellType_boxplots"), h=6)
    ## Composition barplot
    prop_df = memPROP() %>% filter(Tissue == tissue)
    comp_barplot <- plot_composition_bar(prop_df,x_col="Race",sample_col="sample") +
        ggpubr::theme_pubr(base_size = 15) +
        scale_fill_manual(values = cell_colors)
    ## Save image
    save_img(comp_barplot, paste0(dirname, "/composition_barplot"))
}

#### Correlation with qSV ####
for(tissue in c("Caudate", "DLPFC", "Hippocampus", "Dentate Gyrus")){
    dirname = gsub(" ", "_", tolower(tissue))
    qSV_mat <- read.csv(get_qsv_file(tissue), row.names = "X")
    qSV_long <- qSV_mat %>% rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to = "qSV", values_to = "qSV_value")
    ## Bind with qSV table
    est_prop_qsv <- memPROP() %>% filter(Tissue == tissue) %>%
        mutate(RNum = sample) %>% left_join(qSV_long, by = "RNum")
    est_prop_qsv$qSV <- factor(est_prop_qsv$qSV, levels = paste0('PC', 1:13))
    ## Calculate p-values ##
    prop_qSV_fit <- est_prop_qsv %>% group_by(cell_type, qSV) %>%
        do(fitQSV = tidy(lm(prop ~ qSV_value, data = .))) %>%
        unnest(fitQSV) %>% filter(term == "qSV_value") %>%
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf))
    print(prop_qSV_fit %>% count(p.bonf.cat))
    ## Tile plots
    my_breaks <- c(0.05, 0.01, 0.005, 0)
    tile_plot_val <- prop_qSV_fit %>%
        ggplot(aes(x = cell_type, y = qSV, fill = log.p.bonf)) +
        geom_tile(color = "grey") +
        geom_text(aes(label=ifelse(p.bonf.sig,
                                   format(round(log.p.bonf,1), nsmall=1), ""),
                      color=log.p.bonf),
                  size=3, fontface="bold", show.legend=FALSE) +
        scale_color_viridis(option = "magma") +
        scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                           direction=-1) +
        labs(title ="p-values cell-type prop~qSV", x = 'Cell Type',
             color ="p-value Bonf\nsignificance") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ggpubr::theme_pubr(base_size = 15)
    save_img(tile_plot_val, paste0(dirname, "/qSV_prop_fit_tileVal"))
    ## Create scatter plots
    sig_colors2 <- c("#440154","#31688E", "#35B779","black")
    names(sig_colors2) <- levels(prop_qSV_fit$p.bonf.cat)
    est_prop_qsv_fit <- left_join(est_prop_qsv, prop_qSV_fit)
    scatter_plot_cat <- est_prop_qsv_fit %>%
        ggplot(aes(x = qSV_value, y = prop, color = p.bonf.cat))+
        geom_point(size = .4, alpha = .5) +
        facet_grid(cell_type~qSV, scales = "free") +
        scale_color_manual(values = sig_colors2) +
        ggpubr::theme_pubr(base_size = 15)+
        theme(legend.text = element_text(size = 15)) +
        guides(color = guide_legend(override.aes = list(size=5)))
    save_img(scatter_plot_cat, paste0(dirname, "/qSV_cellType_scatter_cat"),
             26, 18)
    scatter_plot_val <- est_prop_qsv_fit %>%
        ggplot(aes(x = qSV_value, y = prop, color = log.p.bonf))+
        geom_point(size = .4, alpha = .6) +
        facet_grid(cell_type~qSV, scales = "free")+
        scale_color_viridis(name = "-log10(p-value Bonf)", option = "magma",
                            direction = -1) +
        ggpubr::theme_pubr(base_size = 15)+
        theme(legend.text = element_text(size = 15)) +
        guides(color = guide_legend(override.aes = list(size=5)))
    save_img(scatter_plot_val, paste0(dirname, "/qSV_cellType_scatter_val"),
             26, 18)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
