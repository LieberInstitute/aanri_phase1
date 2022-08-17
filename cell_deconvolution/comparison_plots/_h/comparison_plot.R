## Cell deconvolution comparison and quality control

library(tidyverse)
library(ggpubr)

get_pheno <- function(races){
    baseloc = "../../../input"
    ancestry <- paste0(baseloc,"/ancestry_structure/structure.out_ancestry",
                       "_proportion_raceDemo_compare")
    fname = paste0(baseloc,"/phenotypes/merged/_m/merged_phenotypes.csv")
    df = data.table::fread(fname) %>% select(-V1) %>%
        filter(Dx %in% c("Control"), Age > 17,
               Race %in% races) %>%
        inner_join(data.table::fread(ancestry),
                   by=c("BrNum"="id", "Race"="group")) %>%
        mutate(Race = gsub("CAUC", "EA", Race))
    return(df)
}
memPHENO <- memoise::memoise(get_pheno)

save_img <- function(image, fn, w, h){
    for(ext in c(".svg", ".pdf")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

prep_data <- function(races){
                                        # Prepare data
    load("../../_m/est_prop_Bisque.v2.Rdata", verbose = TRUE)

    cc = est_prop_bisque$caudate$Est.prop.long %>% 
        inner_join(memPHENO(races), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="Caudate")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>% 
        inner_join(memPHENO(races), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="DLPFC")
    hh = est_prop_bisque$hippo$Est.prop.long %>% 
        inner_join(memPHENO(races), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="Hippocampus")
    gg = est_prop_bisque$dg$Est.prop.long %>% 
        separate(sample, c("sample", "batch")) %>% 
        inner_join(memPHENO(races), by=c("sample"="RNum")) %>%
        mutate_if(is.character, as.factor) %>%
        rename("Proportion"="prop") %>% mutate(Tissue="Dentate Gyrus")
    return(bind_rows(cc, dd, hh, gg))
}
memDF <- memoise::memoise(prep_data)

generate_plots <- function(races, label){
    df <- memDF(races)
                                        # Cell type propotion plots
    bxp = df %>% ggboxplot(x="cell_type", y="Proportion", color="Afr", facet.by="Tissue",
                           panel.labs.font=list(face='bold', size = 14), #palette="npg", 
                           outlier.shape=NA, ylab='Cell Type Proportion', add='jitter', 
                           add.params=list(alpha=0.5), ylim=c(0, 1), xlab="Cell Types", 
                           legend="bottom") +
        font("xy.text", size=12) + font("xy.title", size=16, face="bold") + 
        rotate_x_text(45)
    save_img(bxp, paste0("boxplot_celltypes_ancestry_", label), w=7, h=6)

    bxp = df %>%
        ggscatter(x="Afr", y="Proportion", color="Afr", facet.by=c("cell_type","Region"),
                  panel.labs.font=list(face='bold', size = 14), #palette="npg", 
                  ylab='Cell Type Proportion', add='reg.line', ncol=4, 
                  add.params=list(color="blue", fill="lightgray"), conf.int=TRUE, 
                  cor.coef=TRUE, xlab="Genetic Ancestry") +
        font("xy.text", size=12) + font("xy.title", size=16, face="bold")
    fn = paste0("scatterplot_ancestryBYcelltype_",label)
    save_img(bxp, fn, w=9, h=12)
}

#### MAIN ####
generate_plots(c("AA", "CAUC"), "all")
generate_plots(c("AA"), "AAonly")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
