## Differential expression with limma-voom pipeline
suppressMessages({
    library(SummarizedExperiment)
    library(tidyverse)
})

merge_rse_metrics <- function(rse) {
                                        # Function from jaffelab github
    stopifnot(is(rse, 'RangedSummarizedExperiment'))

    rse$overallMapRate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$overallMapRate, rse$numReads)
    rse$mitoRate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$mitoRate, rse$numMapped)
    rse$rRNA_rate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$rRNA_rate, rse$numMapped)
    rse$totalAssignedGene = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$totalAssignedGene, rse$numMapped)

    rse$numMapped = sapply(rse$numMapped, sum)
    rse$numReads = sapply(rse$numReads, sum)
    rse$numUnmapped = sapply(rse$numUnmapped, sum)
    rse$mitoMapped = sapply(rse$mitoMapped, sum)
    rse$totalMapped = sapply(rse$totalMapped, sum)
    return(rse)
}

save_volcanoPlot <- function(top, label, feature){
    pdf(file=paste0(feature, "/volcanoPlot_", label, ".pdf"), 8, 6)
    with(top, plot(logFC, -log10(P.Value), pch=20, cex=0.6))
    with(subset(top, adj.P.Val<=0.05), points(logFC, -log10(P.Value),
                                              pch=20, col='red', cex=0.6))
    with(subset(top, abs(logFC)>0.50), points(logFC, -log10(P.Value),
                                              pch=20, col='orange', cex=0.6))
    with(subset(top, adj.P.Val<=0.05 & abs(logFC)>0.50),
         points(logFC, -log10(P.Value), pch=20, col='green', cex=0.6))
    dev.off()
}

save_MAplot <- function(top, label, feature){
    pdf(file=paste0(feature, "/MAplot_", label, ".pdf"), 8, 6)
    with(top, plot(AveExpr, logFC, pch=20, cex=0.5))
    with(subset(top, adj.P.Val<0.05),
         points(AveExpr, logFC, col="red", pch=20, cex=0.5))
    dev.off()
}

extract_de <- function(contrast, label, efit, feature){
    top <- limma::topTable(efit, coef=contrast, number=Inf, sort.by="P")
    top$SE <- sqrt(efit$s2.post) * efit$stdev.unscaled
    top.fdr <- top %>% filter(adj.P.Val<=0.05)
    print(paste("Comparison for:", label))
    print(paste('There are:', dim(top.fdr)[1], 'DE features!'))
    data.table::fwrite(top, 
                       file=paste0(feature, "/diffExpr_", label, "_full.txt"), 
                       sep='\t', row.names=TRUE)
    data.table::fwrite(top.fdr, 
                       file=paste0(feature, "/diffExpr_", label, "_FDR05.txt"), 
                       sep='\t', row.names=TRUE)
    save_volcanoPlot(top, label, feature)
    save_MAplot(top, label, feature)
}

prep_data <- function(feature){
    ancestry = "../../../../input/ancestry_structure/structure.out_ancestry_proportion_raceDemo_compare"
    counts_lt = list("genes"="../../../../input/counts/_m/astellas_dg_hg38_rseGene_n263.rda", 
                     "transcripts"="../../../../input/counts/_m/astellas_dg_hg38_rseTx_n263.rda",
                     "exons"="../../../../input/counts/_m/astellas_dg_hg38_rseExon_n263.rda",
                     "junctions"="../../../../input/counts/_m/astellas_dg_hg38_rseJxn_n263.rda")
    tx_file = "../../../../input/counts/_m/transcripts_counts/dg_counts.txt"
    load(counts_lt[[feature]], verbose=TRUE)
    if(exists("rse_gene")){
        rse_df = rse_gene
    } else if (exists("rse_tx")){
        rnames <- sapply(strsplit(rownames(colData(rse_tx)), "_"), "[", 1)
        counts <- data.table::fread(tx_file) %>% column_to_rownames("transcript_id") 
        colnames(counts) <- sapply(strsplit(colnames(counts), "_"), "[", 1)
        counts <- counts %>% select(any_of(rnames)) %>% as.matrix
        annot <- rowData(rse_tx) %>% as.data.frame %>% 
            filter(transcript_id %in% rownames(counts)) %>% as.matrix
        keepIdx <- which(rnames %in% colnames(counts))
        rse_tx  <- rse_tx[, keepIdx]
        rse_df <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                       rowData=annot, colData=colData(rse_tx)) %>%
            as("RangedSummarizedExperiment")
        colnames(rse_df) <- rnames[keepIdx]
    } else if (exists("rse_exon")){
        rse_df = rse_exon
    } else {
        rse_df = rse_jxn
    }
    keepIndex = which((rse_df$Dx %in% c("Control")) & 
                      rse_df$Age > 17 & rse_df$Race %in% c("AA", "CAUC"))
    rse_df = rse_df[, keepIndex]
    rse_df$Dx = factor(rse_df$Dx, levels = c("Control", "Schizo"))
    rse_df$Sex <- factor(rse_df$Sex)
    rse_df <- merge_rse_metrics(rse_df)
    colData(rse_df)$RIN = sapply(colData(rse_df)$RIN,"[",1)
    rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"), "[", 1)
    pheno = colData(rse_df) %>% as.data.frame %>% 
        inner_join(data.table::fread(ancestry), by=c("BrNum"="id", "Race"="group")) %>%
        filter(Afr >= 0.8 | Eur > 0.99)
                                        # Generate DGE list
    x <- edgeR::DGEList(counts=assays(rse_df)$counts[, pheno$RNum], 
                        genes=rowData(rse_df), samples=pheno)
                                        # Filter by expression
    design0 <- model.matrix(~Race, data=x$samples)
    keep.x <- edgeR::filterByExpr(x, design=design0)
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
                                        # Normalize library size
    x <- edgeR::calcNormFactors(x, method="TMM")
    return(x)
}
memPREP <- memoise::memoise(prep_data)

qSV_model <- function(feature){
    x <- memPREP(feature)
    # Design matrix
    mod = model.matrix(~Race + Sex + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + overallMapRate, data = x$samples)
    #colnames(mod) <- gsub("Dx", "", colnames(mod))
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("RaceCAUC", "EA", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Load qSV
    qsv_file = "../../../../input/phenotypes/_m/qSV_dg.csv"
    modQsva <- mod %>% as.data.frame %>% rownames_to_column() %>%
        inner_join(data.table::fread(qsv_file), by=c("rowname"="V1")) %>% 
        rename_all(list(~str_replace_all(., 'PC', 'qPC'))) %>% 
        column_to_rownames("rowname") %>% as.matrix
    return(modQsva)
}
memQSV <- memoise::memoise(qSV_model)

get_voom <- function(feature){
    ### Preform voom
    x <- memPREP(feature)
    modQsva <- memQSV(feature)
    v <- limma::voom(x[, rownames(modQsva)], modQsva)
    return(v)
}
memVOOM <- memoise::memoise(get_voom)

cal_res <- function(feature){
    ### Calculate residuals
    v <- memVOOM(feature)
    null_model <- v$design %>% as.data.frame %>% select(-c("EA")) %>% as.matrix
    fit_res <- limma::lmFit(v, design=null_model)
    res = v$E - ( fit_res$coefficients %*% t(null_model) )
    res_sd = apply(res, 1, sd)
    res_mean = apply(res, 1, mean)
    res_norm = (res - res_mean) / res_sd
    write.table(res_norm, file=paste0(feature, '/residualized_expression.tsv'), 
                sep="\t", quote=FALSE)
}
memRES <- memoise::memoise(cal_res)

fit_voom <- function(feature){
    v       <- memVOOM(feature)
    modQsva <- memQSV(feature)
    fit0    <- limma::lmFit(v, modQsva)
    contr.matrix <- limma::makeContrasts(EAvsAA=EA,
                                         levels=colnames(modQsva))
    fit <- limma::contrasts.fit(fit0, contrasts=contr.matrix)
    esv <- limma::eBayes(fit)
    return(esv)
}
memFIT <- memoise::memoise(fit_voom)

#### MAIN analysis
for(feature in c('genes', 'transcripts', 'junctions', 'exons')){
    dir.create(feature)
    # Preform voom
    v <- memVOOM(feature)
    save(v, file=paste0(feature,'/voomSVA.RData'))
    # Fit model and apply eBayes
    efit = memFIT(feature)
    # Save differential expression
    extract_de(1, "EAvsAA", efit, feature)
    # Calculate residuals
    memRES(feature)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
