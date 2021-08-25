library(tidyverse)

## Functions
save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".pdf", ".png")){
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
    load("../../_m/est_prop_Bisque.v2.Rdata")
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

## Check correlations with cell type proportions
datalist = list()
for(var in c("RIN", "Age", "mitoRate", "rRNA_rate", "overallMapRate")){
    formula = paste0("prop~",var)
    df_fit = memPROP() %>% group_by(Region, cell_type) %>%
        do(lmfit=broom::tidy(lm(formula, data=.))) %>%
        unnest(lmfit) %>% filter(term != "(Intercept)") %>%
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf))
    datalist[[var]] <- df_fit
    ##print(df_fit)
    ##print(df_fit %>% count(p.bonf.cat))
}
bigdf = bind_rows(datalist)
bigdf %>% data.table::fwrite("celltype_proportion_confounders.tsv", sep='\t')

## Tile plot (heatmap)
my_breaks <- c(0.05, 0.01, 0.005, 0)
tile_plot <- bigdf %>%
    ggplot(aes(x = cell_type, y = term, fill = log.p.bonf,
               label = ifelse(p.bonf.sig,
                              format(round(log.p.bonf,1), nsmall=1), ""))) +
    geom_tile(color = "grey") + ggfittext::geom_fit_text(contrast = TRUE) +
    facet_wrap("~Region", nrow=2) +
    viridis::scale_color_viridis(option = "magma") +
    viridis::scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                                direction=-1) +
    labs(x="Cell Type", color="p-value Bonf\nsignificance", y="Covariates") +
    ggpubr::theme_pubr(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_text(face="bold"))

save_img(tile_plot, "tilePlot_confounders_vs_celltypes")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
