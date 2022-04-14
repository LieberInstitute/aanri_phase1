## Calculate residualized expression
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(variancePartition)
    library(BiocParallel)
    library(tidyverse)
    library(glue)
})

set_parallel <- function(){
    parallel_param = 32
    param = SnowParam(parallel_param, "SOCK", progressbar=FALSE)
    register(param)
}

get_cell_prop <- function(){
    celltype_file <- paste0("../../../../cell_deconvolution/_m/",
                            "est_prop_Bisque.v2.Rdata")
    load(celltype_file)
    cc = est_prop_bisque$caudate$Est.prop.long %>%
        mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    dd = est_prop_bisque$dlpfc$Est.prop.long %>%
        mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    hh = est_prop_bisque$hippo$Est.prop.long %>%
        mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    gg = est_prop_bisque$dg$Est.prop.long %>%
        separate(sample, c("sample", "batch")) %>%
        select(-batch) %>% mutate_if(is.character, as.factor) %>%
        rename("proportion"="prop", "RNum"="sample")
    return(bind_rows(cc, dd, hh, gg) %>%
           pivot_wider(names_from="cell_type", values_from="proportion"))
}
memPROP <- memoise::memoise(get_cell_prop)

merge_rse_metrics <- function(rse) {
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

load_counts <- function(region){
    if(region == 'caudate'){
        counts_lt = paste0("../../../../input/counts/_m/caudate_",
                           "brainseq_phase3_hg38_rseGene_merged_n464.rda")
    } else if(region == 'dlpfc'){
        counts_lt = paste0("../../../../input/counts/_m/dlpfc_ribozero_",
                           "brainseq_phase2_hg38_rseGene_merged_n453.rda")
    } else if(region == 'hippocampus'){
        counts_lt = paste0("../../../../input/counts/_m/hippo_",
                           "brainseq_phase2_hg38_rseGene_merged_n447.rda")
    } else {
        counts_lt = "../../../../input/counts/_m/astellas_dg_hg38_rseGene_n263.rda"
    }
    return(list("count_file"=counts_lt))
}

check_dup <- function(df){
    sample <- df %>% select_if(is.numeric)
    variables <- names(sample)
    return(cytominer::correlation_threshold(variables, sample, cutoff=0.95))
}

prep_data <- function(){
    ancestry = paste0("../../../../input/ancestry_structure/structure.",
                      "out_ancestry_proportion_raceDemo_compare")
    datalist1 = list(); datalist2 = list()
    for(region in c('dlpfc', 'caudate', 'hippocampus', 'gyrus')){
        print(region)
        flush.console()
        lt = load_counts(region)
        load(lt[['count_file']])
        rse_df = rse_gene
        keepIndex = which((rse_df$Dx %in% c("Control")) & rse_df$Age > 17 &
                          rse_df$Race %in% c("AA"))
        rse_df = rse_df[, keepIndex]
        rse_df$Sex <- factor(rse_df$Sex)
        rse_df <- merge_rse_metrics(rse_df)
        colData(rse_df)$RIN = sapply(colData(rse_df)$RIN,"[",1)
        rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"), "[", 1)
        pheno = colData(rse_df) %>% as.data.frame %>%
            inner_join(data.table::fread(ancestry), by=c("BrNum"="id", "Race"="group"))
        datalist1[[region]] <- assays(rse_df)$counts[, pheno$RNum] %>% t %>%
                as.data.frame
        datalist2[[region]] <- pheno[, c("Eur", 'Sex', 'Age', 'mitoRate',
                                         'rRNA_rate', 'totalAssignedGene',
                                         'overallMapRate', 'Region', 'BrNum', "RNum")]
    }
    samples <- bind_rows(datalist2) %>% as.data.frame %>%
        inner_join(memPROP(), by="RNum") %>%
        mutate_if(is.numeric, scales::rescale)
    if(length(check_dup(samples)) != 0){
        samples <- samples %>% select(-check_dup(samples))
    }
    subjects = samples$BrNum
    genes = rowData(rse_df) %>% as.data.frame
    counts <- bind_rows(datalist1) %>% t %>% as.data.frame
    colnames(counts) <- samples$RNum
    if(samples %>% has_rownames){
        samples = samples %>% select(-RNum)
    } else {
        samples = samples %>% column_to_rownames("RNum")
    }
    # Generate DGE list
    x <- edgeR::DGEList(counts=counts, genes=genes, samples=samples)
    # Filter by expression
    design0 <- model.matrix(~Region, data=x$samples)
    keep.x <- edgeR::filterByExpr(x, design=design0)
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    # Normalize library size
    x <- edgeR::calcNormFactors(x, method="TMM")
    return(x)
}
memPREP <- memoise::memoise(prep_data)

get_model <- function(){
    covs = c('Sex', 'Age', 'mitoRate', 'rRNA_rate', 'totalAssignedGene',
             'overallMapRate', 'Astro', "Macro", "Micro", "Mural",
             "Oligo", "OPC", "Excit", "Inhib", 'Region', 'BrNum')
    adjust_covs = covs[2:14]
    random_effect = covs[16]
    formula <- glue("~ Eur * Region + ",
                    glue_collapse(adjust_covs, sep = " + ")) %>%
        glue(" + (1|Sex) + (1|{random_effect})")
    return(formula)
}

get_null <- function(){
    covs = c('Sex', 'Age', 'mitoRate', 'rRNA_rate', 'totalAssignedGene',
             'overallMapRate', 'Astro', "Macro", "Micro", "Mural",
             "Oligo", "OPC", "Excit", "Inhib", 'Region', 'BrNum')
    adjust_covs = covs[2:14]
    random_effect = covs[16]
    null_form <- glue("~ ", glue_collapse(adjust_covs, sep = " + ")) %>%
        glue(" + (1|Sex) + (1|{random_effect})")
    return(null_form)
}

get_voom <- function(){
    x <- memPREP()
    vobjDream = voomWithDreamWeights(x, get_model(), x$samples)
    return(vobjDream)
}
memVOOM <- memoise::memoise(get_voom)

cal_res <- function(){
### Calculate residuals
    vobj <- memVOOM()
    residList <- fitVarPartModel(vobj, get_null(), vobj$targets, fxn=residuals)
    do.call(rbind, residList) %>% as.data.frame %>% rownames_to_column("Geneid") %>%
        data.table::fwrite(file='residualized_expression.tsv', sep='\t')
}

#### Main
print("Set parallel!")
set_parallel()
print("Generated residualized expression!")
cal_res()

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
