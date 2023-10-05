## Differential expression with limma-voom pipeline
suppressMessages({
    library(SummarizedExperiment)
    library(tidyverse)
    library(optparse)
    library(here)
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

extract_de <- function(contrast, efit, feature, perm_num){
    top    <- limma::topTable(efit, coef=contrast, number=Inf, sort.by="P")
    top$SE <- sqrt(efit$s2.post) * efit$stdev.unscaled
    fname  <- paste0(feature, "/de_model_ancestry_perm_", perm_num, ".txt")
    top %>%
        data.table::fwrite(file=fname, sep='\t', row.names=TRUE)
}

prep_data <- function(feature, nsize){
    ancestry = here("input/ancestry_structure",
                    "structure.out_ancestry_proportion_raceDemo_compare")
    counts_lt = list("genes"=here("input/counts/_m",
                                  "hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda"), 
                     "transcripts"=here("input/counts/_m",
                                        "hippo_brainseq_phase2_hg38_rseTx_merged_n447.rda"),
                     "exons"=here("input/counts/_m",
                                  "hippo_brainseq_phase2_hg38_rseExon_merged_n447.rda"),
                     "junctions"=here("input/counts/_m",
                                      "hippo_brainseq_phase2_hg38_rseJxn_merged_n447.rda"))
    tx_file = here("input/counts/_m/transcripts_counts/hippo_counts.txt")
    load(counts_lt[[feature]], verbose=TRUE)
    if(exists("rse_gene")){
        rse_df = rse_gene
    } else if (exists("rse_tx")){
        counts <- data.table::fread(tx_file) %>% column_to_rownames("transcript_id") %>%
            select(any_of(colnames(rse_tx))) %>% as.matrix
        annot <- rowData(rse_tx) %>% as.data.frame %>% 
            filter(transcript_id %in% rownames(counts)) %>% as.matrix
        keepIdx <- which(colnames(rse_tx) %in% colnames(counts))
        rse_tx  <- rse_tx[, keepIdx]
        rse_df <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                       rowData=annot, colData=colData(rse_tx)) %>%
            as("RangedSummarizedExperiment")
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
        filter(Afr >= 0.8 | Eur > 0.99) %>%
        group_by(Race) %>% sample_n(size=nsize)
                                        # Generate DGE list
    x <- edgeR::DGEList(counts=assays(rse_df)$counts[, pheno$RNum], 
                        genes=rowData(rse_df), samples=pheno)
                                        # Filter by expression
    design0 <- model.matrix(~Race, data=x$samples)
    keep.x <- edgeR::filterByExpr(x, design=design0)
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    x <- x[keep.x, , keep.lib.sizes=FALSE]
                                        # Normalize library size
    x <- edgeR::calcNormFactors(x, method="TMM")
    return(x)
}
memPREP <- memoise::memoise(prep_data)

qSV_model <- function(feature, nsize){
    x <- memPREP(feature, nsize)
    # Design matrix
    mod = model.matrix(~Race + Sex + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + overallMapRate, data = x$samples)
    #colnames(mod) <- gsub("Dx", "", colnames(mod))
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("RaceCAUC", "EA", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Load qSV
    qsv_file = here("input/phenotypes/_m/qSV_hippo.csv")
    modQsva <- mod %>% as.data.frame %>% rownames_to_column() %>%
        inner_join(data.table::fread(qsv_file), by=c("rowname"="V1")) %>% 
        rename_all(list(~str_replace_all(., 'PC', 'qPC'))) %>% 
        column_to_rownames("rowname") %>% as.matrix
    return(modQsva)
}
memQSV <- memoise::memoise(qSV_model)

get_voom <- function(feature, nsize){
    ### Preform voom
    x <- memPREP(feature, nsize)
    modQsva <- memQSV(feature, nsize)
    v <- limma::voom(x[, rownames(modQsva)], modQsva)
    return(v)
}
memVOOM <- memoise::memoise(get_voom)

fit_voom <- function(feature, nsize){
    v            <- memVOOM(feature, nsize)
    modQsva      <- memQSV(feature, nsize)
    fit0         <- limma::lmFit(v, modQsva)
    contr.matrix <- limma::makeContrasts(EAvsAA=EA,
                                         levels=colnames(modQsva))
    fit <- limma::contrasts.fit(fit0, contrasts=contr.matrix)
    esv <- limma::eBayes(fit)
    return(esv)
}
memFIT <- memoise::memoise(fit_voom)

#### MAIN analysis
                                        # Parser
option_list <- list(
    make_option(c("-p", "--perm_num"), type="integer", default=1,
                help="permutation number [default=%default]"),
    make_option(c("-n", "--n_size"), type="integer", default=61,
                help="stratified samples size number [default=%default]")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)

perm <- opt$perm_num; nsize <- opt$n_size
                                        # Permutate target phenotype
perm_num <- str_pad(perm, 2, pad="0")
seed     <- perm + 13131313
                                        # Set permutation seed
set.seed(seed)

for(feature in c('genes', 'transcripts', 'junctions', 'exons')){
    dir.create(feature)
                                        # Preform voom
    v    <- memVOOM(feature, nsize)
    save(v, file=paste0(feature,'/voomSVA_perm_',perm_num,'.RData'))
                                        # Fit model and apply eBayes
    efit <- memFIT(feature, nsize)
                                        # Save differential expression
    extract_de(1, efit, feature, perm_num)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
