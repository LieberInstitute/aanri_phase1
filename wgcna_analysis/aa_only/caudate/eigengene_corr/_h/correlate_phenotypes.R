## Exploring the eign values and correlate with genetic ancestry
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(here)
    library(dplyr)
    library(ggpubr)
})

get_pheno <- function(){
    fn <- here("input/phenotypes/merged/_m/merged_phenotypes.csv")
    return(data.table::fread(fn))
}

get_ancestry <- function(){
    fn <- here("input/ancestry_structure/",
               "structure.out_ancestry_proportion_raceDemo_compare")
    return(data.table::fread(fn))
}

load_eigenvalues <- function(){
    return(data.table::fread("../../_m/eigengenes.csv"))
}
memEIGEN <- memoise::memoise(load_eigenvalues)

get_modules <- function(){
    return(memEIGEN() %>% select(-V1) %>% colnames)
}

merge_data <- function(){
    dt <- memEIGEN() %>% left_join(get_pheno(), by="V1") %>%
        inner_join(get_ancestry(), by=c("BrNum"="id"))
    return(dt)
}
memDT <- memoise::memoise(merge_data)

main <- function(){
                                        # Linear model
    pvals <- c(); beta <- c()
    for(mod in get_modules()){
        model  <- paste0("Eur ~ ", mod)
        fitted <- lm(model, data=memDT())
        res    <- summary(fitted)
        pvals  <- c(pvals, res$coefficients[mod, "Pr(>|t|)"])
        beta   <- c(beta, res$coefficients[mod, "Estimate"])
    }
    fdr <- p.adjust(pvals, method="fdr")
    df1 <- data.frame("Modules"=get_modules(), "Beta"=beta,
                      "Pvalue"=pvals, "FDR"=fdr)
    print(df1 %>% filter(Pvalue < 0.05))
    df1 %>% mutate(Tissue="Caudate") %>%
        data.table::fwrite("eigen_correlation_ancestry.tsv", sep='\t')
                                        # Pearson correlation
    pvals <- c(); est <- c()
    for(mod in get_modules()){
        res <- cor.test(memDT()[["Eur"]], memDT()[[mod]], method="pearson")
        pvals = c(pvals, res$p.value)
        est = c(est, res$estimate[[1]])
    }
    fdr <- p.adjust(pvals, method="fdr")
    df2 = data.frame("Modules"=get_modules(), "Rho"=est, "Pvalue"=pvals, "FDR"=fdr)
    print(df2 %>% filter(Pvalue < 0.05))
}

#### MAIN
main()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
