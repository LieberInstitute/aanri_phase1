## Differential expression with limma-voom pipeline
suppressMessages({
    library(tidyverse)
})

extract_de <- function(contrast, label, efit, feature){
    top <- limma::topTable(efit, coef=contrast, number=Inf, sort.by="P")
    top$SE <- sqrt(efit$s2.post) * efit$stdev.unscaled
    top.fdr <- top %>% filter(adj.P.Val<=0.05)
    print(paste("Comparison for:", label))
    print(paste('There are:', dim(top.fdr)[1], 'DE features!'))
    top %>% rownames_to_column("phenotype_id") %>%
        data.table::fwrite(file=paste0(feature, "/diffExpr_", label, "_full.txt"), 
                           sep='\t', row.names=FALSE)
    top.fdr %>% rownames_to_column("phenotype_id") %>%
        data.table::fwrite(file=paste0(feature, "/diffExpr_", label, "_FDR05.txt"), 
                           sep='\t', row.names=FALSE)
}

load_model <- function(feature){
    fn <- here::here("differential_analysis/hippocampus/_m/",
                     feature,"/voomSVA.RData")
    load(fn)
    design0 <- v$targets %>% select(BrNum, RNum)
    return(v$design %>% as.data.frame %>% rownames_to_column("RNum") %>%
        inner_join(design0, by="RNum") %>% column_to_rownames("BrNum") %>%
        select(Intercept, EA) %>% as.matrix)
}

fit_model <- function(feature){
    # Load data and design matrix
    res_file <- paste0("../../_m/",feature,".predicted_expression.hippocampus.txt.gz")
    res_expr <- data.table::fread(res_file) %>% column_to_rownames("phenotype_id")
    mod      <- load_model(feature)
    # Fit model
    sample_ids   <- intersect(rownames(mod), colnames(res_expr))
    logCPM       <- res_expr %>% select(any_of(sample_ids))
    design       <- mod[sample_ids, ]
    fit0         <- limma::lmFit(logCPM, design)
    contr.matrix <- limma::makeContrasts(EAvsAA=EA,levels=colnames(design))
    fit <- limma::contrasts.fit(fit0, contrasts=contr.matrix)
    esv <- limma::eBayes(fit)
    return(esv)
}
memFIT <- memoise::memoise(fit_model)

#### MAIN analysis
for(feature in c('genes', 'transcripts', 'junctions', 'exons')){
    dir.create(feature)
    # Fit model and apply eBayes
    efit = memFIT(feature)
    # Save differential expression
    extract_de(1, "EAvsAA", efit, feature)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
