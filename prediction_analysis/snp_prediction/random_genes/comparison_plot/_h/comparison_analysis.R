## Comparison of models and plot

library(tidyverse)
library(ggpubr)

save_plot <- function(p, fn, w=7, h=6){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}


get_ml_summary <- function(fn){
    ml_df = data.table::fread(fn) %>% mutate_at("fold", as.character) %>%
        select(fold, n_features, n_redundant, starts_with("test_score_r2")) %>%
        pivot_longer(-fold) %>% group_by(name) %>%
        summarise(Mean=mean(value), Median=median(value), Std=sd(value), .groups = "keep")
    return(ml_df)
}


get_metrics <- function(filename, tissue){
    datalist = list()
    for(fn in Sys.glob(filename)){
        gene_id = str_extract(fn, "ENSG\\d+_\\d+")
        dat <- get_ml_summary(fn)
        dat["Geneid"] = gene_id
        datalist[[gene_id]] <- dat
    }
    ml_df <- bind_rows(datalist)
    ml_df["Tissue"] = tissue
    return(ml_df)
}


generate_raffe_metrics <- function(){
    filename1 = "../../raffe/_m/caudate/*/raffe_10folds.txt"
    filename2 = "../../raffe/_m/dentateGyrus/*/raffe_10folds.txt"
    filename3 = "../../raffe/_m/dlpfc/*/raffe_10folds.txt"
    filename4 = "../../raffe/_m/hippocampus/*/raffe_10folds.txt"

    raffe_cc <- get_metrics(filename1, "Caudate")
    raffe_gg <- get_metrics(filename2, "Dentate Gyrus")
    raffe_dd <- get_metrics(filename3, "DLPFC")
    raffe_hh <- get_metrics(filename4, "Hippocampus")

    raffe = bind_rows(raffe_cc, raffe_gg, raffe_dd, raffe_hh) %>%
        as.data.frame %>% mutate_if(is.character, as.factor) %>%
        filter(name == "test_score_r2") %>%
        mutate("Model" = "Random Forest") %>% select(-"name")
    ##print(dim(raffe))
    ##print(raffe %>% head(2))
    return(raffe)
}


generate_enet_metrics <- function(){
    filename1 = "../../enet/_m/caudate/*/enet_rfe_10folds.txt"
    filename2 = "../../enet/_m/dentateGyrus/*/enet_rfe_10folds.txt"
    filename3 = "../../enet/_m/dlpfc/*/enet_rfe_10folds.txt"
    filename4 = "../../enet/_m/hippocampus/*/enet_rfe_10folds.txt"

    enet_cc <- get_metrics(filename1, "Caudate")
    enet_gg <- get_metrics(filename2, "Dentate Gyrus")
    enet_dd <- get_metrics(filename3, "DLPFC")
    enet_hh <- get_metrics(filename4, "Hippocampus")

    enet = bind_rows(enet_cc, enet_gg, enet_dd, enet_hh) %>%
        as.data.frame %>% mutate_if(is.character, as.factor) %>%
        filter(name == "test_score_r2") %>%
        mutate("Model" = "Elastic Net") %>% select(-"name")
    ##dim(enet)
    ##enet %>% head(2)
    return(enet)
}


generate_partialR2_metrics <- function(){
    filename1 = "../../partial_r2/_m/caudate/raffe_partial_r2.tsv"
    filename2 = "../../partial_r2/_m/dentateGyrus/raffe_partial_r2.tsv"
    filename3 = "../../partial_r2/_m/dlpfc/raffe_partial_r2.tsv"
    filename4 = "../../partial_r2/_m/hippocampus/raffe_partial_r2.tsv"
    
    partial_cc <- data.table::fread(filename1) %>%
        mutate(Tissue="Caudate")
    partial_gg <- data.table::fread(filename2) %>%
        mutate(Tissue="Dentate Gyrus")
    partial_dd <- data.table::fread(filename3) %>%
        mutate(Tissue="DLPFC")
    partial_hh <- data.table::fread(filename4) %>%
        mutate(Tissue="Hippocampus")

    partial <- bind_rows(partial_cc, partial_gg, partial_dd, partial_hh) %>%
        as.data.frame %>% mutate_if(is.character, as.factor) %>%
        mutate("Model" = "Partial R2")
    return(partial)
}


main <- function(){
    # Partial R2 plots
    bxp <- generate_partialR2_metrics() %>%
        ggboxplot(x="Tissue", y="Partial_R2", fill="Tissue", add="jitter",
                  palette="npg", xlab="", ylab="Partial R2", legend="None",
                  add.params = list(alpha=0.5), ylim=c(-0.25, 1.25),
                  panel.labs.font=list(face="bold", size=14)) +
        rotate_x_text(45) + font("xy.title", size=18, face="bold") +
        font("xy.text", size=16) + font("legend.text", size=16)
    save_plot(bxp, "summary_boxplots_partial_r2", 5, 5)
    # Recursive elimination
    bxp = bind_rows(generate_raffe_metrics(), generate_enet_metrics()) %>%
        ggboxplot(x="Tissue", y="Median", fill="Tissue", add="jitter", 
                  facet.by="Model", palette="npg", ylim=c(-0.75, 1.25), 
                  ylab="Median R2", xlab="", legend="None",
                  add.params=list(alpha=0.5),
                  panel.labs.font=list(face='bold', size = 14)) + 
        rotate_x_text(45) + font("xy.title", size=18, face="bold") + 
        font("xy.text", size=16) + font("legend.text", size=16)
    save_plot(bxp, "summary_boxplots_r2_2methods", 6, 5)
}

options(bitmapType='cairo')
main()
