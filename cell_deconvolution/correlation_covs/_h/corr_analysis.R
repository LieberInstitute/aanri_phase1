library(tidyverse)

## Functions
save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".pdf", ".svg")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_pheno <- function(region){
    baseloc <- paste0("/dcs04/lieber/statsgen/jbenjami/projects/",
                      "aanri_phase1/differential_analysis/")
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
        mutate(across(where(is.factor), as.numeric)) %>%
        mutate(across(where(is.logical), as.numeric)) %>%
        select(Afr,Sex,Age,RIN,mitoRate,rRNA_rate,totalAssignedGene,
               overallMapRate, concordMapRate)
    rownames(dfy) <- rownames(v$targets)
    dfy = dfy %>% rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
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

organize_variables <- function(){
    df0 = memPROP()
    dfx = df0 %>% filter(!grepl("qPC", Covariate)) %>%
        mutate(Covariate=droplevels(Covariate))
    dfy = df0 %>% filter(grepl("qPC", Covariate)) %>%
        mutate(Covariate=droplevels(Covariate))
    dfy$Covariate = factor(dfy$Covariate, levels = paste0("qPC", 1:20))
    levels(dfy$Covariate) <- paste0("qSV", 1:20)
    return(bind_rows(dfx,dfy))
}
memDT <- memoise::memoise(organize_variables)

fit_model <- function(){
    datalist <- list()
    for(tissue in c("Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus")){
        est_fit <- memDT() %>% filter(Tissue == tissue) %>%
            group_by(cell_type, Covariate) %>%
            do(fitEST = broom::tidy(lm(prop ~ Variable, data = .))) %>%
            unnest(fitEST) %>% filter(term != "(Intercept)") %>%
            mutate(p.bonf = p.adjust(p.value, "bonf"),
                   p.bonf.sig = p.bonf < 0.05,
                   p.bonf.cat = cut(p.bonf, breaks=c(1,0.05,0.01,0.005,0),
                                    labels=c("<= 0.005","<= 0.01",
                                             "<= 0.05","> 0.05"),
                                    include.lowest = TRUE),
                   p.fdr = p.adjust(p.value, "fdr"),
                   log.p.bonf = -log10(p.bonf+10**(-300)))
        print(est_fit %>% count(p.bonf.cat))
        datalist[[tissue]] <- est_fit %>% mutate(Tissue = tissue)
    }
    return(bind_rows(datalist))
}
memFIT <- memoise::memoise(fit_model)

tile_plot <- function(){
    my_breaks <- c(0.05, 0.01, 0.005, 0); limits = c(0, 100)    
    tile_plot <- memFIT() %>%
        ggplot(aes(x = cell_type, y = Covariate, fill = log.p.bonf,
                   label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                  nsmall=1), ""))) +
        geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast=TRUE, aes(fontface="bold")) +
        facet_wrap("~Tissue", nrow=2, scales="free") +
        viridis::scale_color_viridis(option="magma") +
        viridis::scale_fill_viridis(name="-log10(p-value Bonf)",
                                    option="magma", direction=-1) +
        labs(x="Cell Type", y="Covariates") +
        ggpubr::theme_pubr(base_size=15) +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              panel.grid=element_blank(),
              axis.title=element_text(size=18, face="bold"),
              strip.text=element_text(face="bold", size=16))
    save_img(tile_plot, "tilePlot_celltypes_covariates", 13, 13)
}

#### MAIN
                                        # Summaize linear model
memFIT() %>%
    data.table::fwrite("celltype_proportion_covariates.tsv", sep='\t')

                                        # Tile plot (heatmap)
tile_plot()

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
