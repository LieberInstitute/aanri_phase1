## Plot explained variance by feature across brain regions
library(dplyr)
library(ggpubr)

tissue_map <- function(tissue){
    return(list("Caudate"="caudate", "Dentate Gyrus"="dentateGyrus",
                "DLPFC"="dlpfc", "Hippocampus"="hippocampus")[[tissue]])
}

get_rsq <- function(feature){
    datalist = list()
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        fn = paste0("../../../",tissue_map(tissue),"/feature_pca/_m/",
                    tolower(feature),"_variance_explained_rsq_1000.csv")
        dfx = data.table::fread(fn) %>% mutate(Tissue=tissue)
        datalist[[tissue_map(tissue)]] <- dfx
    }
    dt = data.table::rbindlist(datalist) %>% mutate(Feature=feature)
    return(dt)
}

get_data <- function(){
    datalist = list()
    for(feature in c("Genes", "Transcripts", "Exons", "Junctions")){
        datalist[[feature]] <- get_rsq(feature)
    }
    return(data.table::rbindlist(datalist))
}

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

plot_scatter <- function(){
    sca <- get_data() %>%
        ggscatter(x="DEGs", y="Rsq", color="Feature", facet.by="Tissue",
                  palette="npg", ylab="Explained Variance\n(PC1, Rsq)",
                  xlab="Ancestry-associated DE Features",
                  ggtheme=theme_pubr(base_size=20, border=F)) +
        theme(axis.title=element_text(face="bold"),
              axis.text=element_text(face="bold"))
    save_plot(sca, "explained_variance", 10, 10)
}

plot_scatter()

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
