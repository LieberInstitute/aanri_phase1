## Author: Louise Huuki
## Edited by KJ Benjamin

suppressPackageStartupMessages({
    library(sessioninfo)
    library(tidyverse)
    library(viridis)
    library(broom)
    library(here)
})

## Functions
get_qsv_file <- function(tissue){
    qsv_lt = list("Caudate"=here("input/phenotypes/_m/qSV_caudate.csv"),
                  "DLPFC"=here("input/phenotypes/_m/qSV_dlpfc.csv"),
                  "Hippocampus"=here("input/phenotypes/_m/qSV_hippo.csv"),
                  "Dentate Gyrus"=here("input/phenotypes/_m/qSV_dg.csv"))
    return(qsv_lt[[tissue]])
}

save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".svg", ".pdf")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

tissue_map <- function(region){
    return(list("Caudate"="caudate", "Dentate Gyrus"="dentateGyrus",
                "DLPFC"="dlpfc", "Hippocampus"="hippocampus")[[region]])
}

get_pheno <- function(region){
    baseloc <- here("differential_analysis/")
    ## Load voom normalized data
    load(paste0(baseloc, region, "/_m/genes/voomSVA.RData"))
    dfx = v$design %>% as.data.frame %>% select(starts_with("qPC")) %>%
        rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    dfx$Covariate <- factor(dfx$Covariate, levels = paste0('qPC', 1:20))
    dfy = v$targets %>% as.data.frame %>%
        mutate(concordMapRate=mapply(function(r, n) {sum(r*n)/sum(n)},
                                     concordMapRate, numMapped)) %>%
        mutate(across(where(is.character), as.factor)) %>%
        select(Afr,Race,Sex,Age,RIN,mitoRate,rRNA_rate,
               totalAssignedGene, overallMapRate, concordMapRate)
    rownames(dfy) <- rownames(v$targets)
    dfy = dfy %>% rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable",
                     values_transform = list(Variable = as.numeric))
    return(bind_rows(dfy, dfx) %>% mutate_if(is.character, as.factor))
}
memPHENO <- memoise::memoise(get_pheno)

get_cell_prop <- function(){
    ## Load Bisque Estimated Props
    load("../../_m/est_prop_Bisque.v2.Rdata")
    cc = est_prop_bisque$caudate$Est.prop.long %>%
        inner_join(memPHENO("caudate"), by=c("sample"="RNum")) %>%
        mutate(Tissue="Caudate")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>%
        inner_join(memPHENO("dlpfc"), by=c("sample"="RNum")) %>%
        mutate(Tissue="DLPFC")
    hh = est_prop_bisque$hippo$Est.prop.long %>%
        inner_join(memPHENO("hippocampus"), by=c("sample"="RNum")) %>%
        mutate(Tissue="Hippocampus")
    gg = est_prop_bisque$dg$Est.prop.long %>%
        separate(sample, c("sample", "batch")) %>%
        inner_join(memPHENO("dentateGyrus"), by=c("sample"="RNum")) %>%
        mutate(Tissue="Dentate Gyrus") %>% select(-batch)
    df = bind_rows(cc, dd, hh, gg)
    return(df)
}
memPROP <- memoise::memoise(get_cell_prop)

#### Load colors and plotting functions ####
cell_types = c("Astro", "Endo", "Micro", "Mural", "Oligo",
               "OPC", "Tcell", "Excit", "Inhib")
cell_colors <- DeconvoBuddies::create_cell_colors(cell_types=cell_types,pallet = "classic")
cell_colors <- cell_colors[c("Astro", "Endo", "Micro", "Mural", "Oligo",
                             "OPC", "Tcell", "Excit", "Inhib")]
#### Comparison plots ####
                                        # Boxplot
group_boxPlot <- memPROP() %>% #filter(Tissue == tissue) %>%
    ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
    geom_boxplot() + labs(x = "Cell Type", y = "Proportion") +
    facet_wrap("~Tissue") + ggpubr::theme_pubr(base_size = 15) +
    scale_fill_manual(values = cell_colors, guide = "none") +
    theme(strip.text = element_text(face="bold"),
          axis.text.x = element_text(angle=45,hjust=1),
          axis.title = element_text(face="bold"))
                                        # Save image
save_img(group_boxPlot, paste0("cellType_boxplots_4regions"), h=9, w=9)

##                                         # Composition barplot
## for(tissue in c("Caudate", "DLPFC", "Hippocampus", "Dentate Gyrus")){
##     dirname = gsub(" ", "_", tolower(tissue))
##     dir.create(dirname)
##     prop_df = memPROP() %>% filter(Tissue == tissue)
##     comp_barplot <- DeconvoBuddies::plot_composition_bar(prop_df,#x_col="Race",
##                                                          sample_col="sample") +
##         ggpubr::theme_pubr(base_size = 15) +
##         scale_fill_manual(values = cell_colors)
##                                         # Save image
##     save_img(comp_barplot, paste0(dirname, "/composition_barplot"))
## }

#### Correlation with qSV ####
for(tissue in c("Caudate", "DLPFC", "Hippocampus", "Dentate Gyrus")){
    dirname = gsub(" ", "_", tolower(tissue))
    qSV_mat <- read.csv(get_qsv_file(tissue), row.names = "X")
    qSV_long <- qSV_mat %>% rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to = "qSV", values_to = "qSV_value")
    ## Bind with qSV table
    est_prop_qsv <- memPROP() %>% filter(Tissue == tissue) %>%
        mutate(RNum = sample) %>% left_join(qSV_long, by = "RNum")
    est_prop_qsv$qSV <- factor(est_prop_qsv$qSV, levels = paste0('PC', 1:15))
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
    tile_plot_val <- prop_qSV_fit %>% mutate(Tissue = tissue) %>%
        ggplot(aes(x = cell_type, y = qSV, fill = log.p.bonf,
                   label=ifelse(p.bonf.sig,
                                format(round(log.p.bonf,1), nsmall=1), ""))) +
        geom_tile(color = "grey") + ggfittext::geom_fit_text(contrast = TRUE) +
        facet_wrap("~Tissue") + scale_color_viridis(option = "magma") +
        scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                           direction=-1, limits=c(0,50)) +
        labs(x = 'Cell Type', color ="p-value Bonf\nsignificance") +
        ggpubr::theme_pubr(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.title = element_text(face="bold", size=18),
              strip.text = element_text(face="bold"))
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
