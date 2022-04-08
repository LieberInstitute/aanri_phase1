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
        fn = paste0("../../../",tissue_map(tissue),"/random_pca/_m/",
                    tolower(feature),
                    "_variance_explained_rsq_100_permutations.csv")
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

plot_boxplot <- function(){
    bxp <- get_data() %>%
        ggboxplot(x="Feature", y="Rsq", color="Feature", facet.by="Tissue",
                  palette="npg", ylab="Explained Variance\n(PC1, Rsq)",
                  xlab="100 Random Features", add="jitter",
                  panel.labs.font=list(face='bold'), outlier.shape=NA,
                  add.params=list(alpha=0.5), ylim=c(0, 0.5),
                  ggtheme=theme_pubr(base_size=20, border=F)) +
        theme(axis.title=element_text(face="bold"),
              axis.text=element_text(face="bold")) +
        rotate_x_text(45)
    save_plot(bxp, "random_explained_variance", 7, 7)
}

plot_boxplot()

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
