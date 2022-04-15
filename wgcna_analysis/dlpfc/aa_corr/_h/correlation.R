## Exploring eigen values and correlation with genetic ancestry for AA only.

suppressPackageStartupMessages({
     library(dplyr)
    library(ggpubr)
})

get_pheno <- function(){
    fn = "../../../../input/phenotypes/merged/_m/merged_phenotypes.csv"
    return(data.table::fread(fn))
}
memPHENO <- memoise::memoise(get_pheno)

get_ancestry <- function(){
    fn = paste0("../../../../input/ancestry_structure/structure.",
                "out_ancestry_proportion_raceDemo_compare")
    return(data.table::fread(fn))
}

load_eigenvalues <- function(){
    return(data.table::fread("../../_m/eigengenes.csv"))
}

merge_data <- function(){
    return(load_eigenvalues() %>% left_join(memPHENO(), by="V1") %>%
           inner_join(get_ancestry(), by=c("BrNum"="id")))
}
memDF <- memoise::memoise(merge_data)

subset_data <- function(){
    return(memDF() %>% filter(Race == "AA"))
}

save_plot <- function(p, fn, w, h){
    for(ext in c('.png', '.pdf')){
        ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

#### MAIN
lf <- logr::log_open("summary.log", logdir=FALSE, autolog=FALSE)
                                        # AA + EA
pvals = c(); est = c(); mods = c()
for(mod in load_eigenvalues() %>% select(-V1) %>% colnames){
    res = cor.test(memDF()[["Eur"]], memDF()[[mod]], method="pearson")
    pvals = c(pvals, res$p.value)
    est = c(est, res$estimate[[1]])
    mods = c(mods, mod)
}
fdr <- p.adjust(pvals, method="fdr")
df1 = data.frame("Modules"=mods, "Rho"=est, "Pvalue"=pvals, "FDR"=fdr)
logr::log_print("Significantly correlated modules: AA+EA")
logr::log_print(df1 %>% filter(Pvalue < 0.05) %>% arrange(Pvalue))
                                        # AA only
pvals = c(); est = c(); mods = c()
for(mod in load_eigenvalues() %>% select(-V1) %>% colnames){
    res = cor.test(subset_data()[["Eur"]], subset_data()[[mod]],
                   method="pearson")
    pvals = c(pvals, res$p.value)
    est = c(est, res$estimate[[1]])
    mods = c(mods, mod)
}
fdr <- p.adjust(pvals, method="fdr")
df2 <- data.frame("Modules"=mods, "Rho"=est, "Pvalue"=pvals, "FDR"=fdr)
logr::log_print("Correlated modules: AA only")
logr::log_print(df2 %>% arrange(Pvalue) %>% head(10))
                                        # Plotting
df  <- data.frame("Modules"=df1$Modules, "All Samples"=df1$Rho,
                  "AA only"=df2$Rho, "logP"=-log2(df1$Pvalue))
df %>% data.table::fwrite("sample_correlation.tsv", sep='\t')

sca <- ggscatter(df, x="All.Samples", y="AA.only", color="logP", size=1,
                 xlab="All Samples (rho)", ylab="AA only (rho)", add="reg.line",
                 panel.labs.font=list(face="bold"), conf.int=TRUE,cor.coef=TRUE,
                 add.params=list(color="blue",fill="lightgray"),cor.coef.size=3,
                 cor.method="pearson", cor.coeff.args=list(label.sep="\n"),
                 ggtheme=ggpubr::theme_pubr(base_size=15))
save_plot(sca, "module_effect_correlation", 6, 6)
#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
