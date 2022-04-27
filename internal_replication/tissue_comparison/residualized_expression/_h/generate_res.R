## Calculate residualized expression
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(variancePartition)
    library(BiocParallel)
    library(tidyverse)
    library(glue)
    library(sva)
})

set_parallel <- function(feature){
    parallel_param = list("genes"=32, "transcripts"=16,
                          "exons"=10, "junctions"=12)
    param = SnowParam(parallel_param[[feature]],
                      "SOCK", progressbar=TRUE)
    register(param)
    return(param)
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

load_counts <- function(region, feature){
    if(region == 'caudate'){
        counts_lt = list("genes"=paste0("../../../../input/counts/_m/",
                                        "caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda"),
                         "transcripts"=paste0("../../../../input/counts/_m/",
                                              "caudate_brainseq_phase3_hg38_rseTx_merged_n464.rda"),
                         "exons"=paste0("../../../../input/counts/_m/",
                                        "caudate_brainseq_phase3_hg38_rseExon_merged_n464.rda"),
                         "junctions"=paste0("../../../../input/counts/_m/",
                                            "caudate_brainseq_phase3_hg38_rseJxn_merged_n464.rda"))
        tx_file = "../../../../input/counts/_m/transcripts_counts/caudate_counts.txt"
    } else if(region == 'dlpfc'){
        counts_lt = list("genes"=paste0("../../../../input/counts/_m/",
                                        "dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n453.rda"),
                         "transcripts"=paste0("../../../../input/counts/_m/",
                                              "dlpfc_ribozero_brainseq_phase2_hg38_rseTx_merged_n453.rda"),
                         "exons"=paste0("../../../../input/counts/_m/",
                                        "dlpfc_ribozero_brainseq_phase2_hg38_rseExon_merged_n453.rda"),
                         "junctions"=paste0("../../../../input/counts/_m/",
                                            "dlpfc_ribozero_brainseq_phase2_hg38_rseJxn_merged_n453.rda"))
        tx_file = "../../../../input/counts/_m/transcripts_counts/dlpfc_counts.txt"
    } else if(region == 'hippocampus'){
        counts_lt = list("genes"=paste0("../../../../input/counts/_m/",
                                        "hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda"),
                         "transcripts"=paste0("../../../../input/counts/_m/",
                                              "hippo_brainseq_phase2_hg38_rseTx_merged_n447.rda"),
                         "exons"=paste0("../../../../input/counts/_m/",
                                        "hippo_brainseq_phase2_hg38_rseExon_merged_n447.rda"),
                         "junctions"=paste0("../../../../input/counts/_m/",
                                            "hippo_brainseq_phase2_hg38_rseJxn_merged_n447.rda"))
        tx_file = "../../../../input/counts/_m/transcripts_counts/hippo_counts.txt"
    } else {
        counts_lt = list("genes"="../../../../input/counts/_m/astellas_dg_hg38_rseGene_n263.rda",
                         "transcripts"="../../../../input/counts/_m/astellas_dg_hg38_rseTx_n263.rda",
                         "exons"="../../../../input/counts/_m/astellas_dg_hg38_rseExon_n263.rda",
                         "junctions"="../../../../input/counts/_m/astellas_dg_hg38_rseJxn_n263.rda")
        tx_file = "../../../../input/counts/_m/transcripts_counts/dg_counts.txt"
    }
    return(list("count_file"=counts_lt[[feature]], "tx_file"=tx_file))
}

get_juncs <- function(region){
    lt = load_counts(region, "junctions")
    load(lt[['count_file']])
    rse_df = rse_jxn
    genes = rowData(rse_df) %>% as.data.frame
    return(rownames(genes))
}

get_overlapping_juncs <- function(){
    common_genes = intersect(intersect(intersect(get_juncs("caudate"),
                                                 get_juncs("dlpfc")),
                                       get_juncs("hippocampus")),
                             get_juncs("gyrus"))
    return(common_genes)
}
memJuncs <- memoise::memoise(get_overlapping_juncs)

check_dup <- function(df){
    sample <- df %>% select_if(is.numeric)
    variables <- names(sample)
    return(cytominer::correlation_threshold(variables, sample, cutoff=0.95))
}

prep_data <- function(feature){
    ancestry = paste0("../../../../input/ancestry_structure/structure.",
                      "out_ancestry_proportion_raceDemo_compare")
    datalist1 = list(); datalist2 = list()
    for(region in c('dlpfc', 'caudate', 'hippocampus', 'gyrus')){
        print(region)
        flush.console()
        lt = load_counts(region, feature)
        load(lt[['count_file']])
        tx_file = lt[['tx_file']]
        if (exists("rse_gene")){
            rse_df = rse_gene
        } else if (exists("rse_tx")){
            if(region == 'caudate' && exists("geneid")){
                counts <- data.table::fread(tx_file) %>%
                    filter(gencodeTx %in% geneid) %>%
                    column_to_rownames("gencodeTx") %>%
                    select(-c(txLength, gencodeID, Symbol, gene_type)) %>%
                    select(colnames(rse_tx)) %>% as.matrix
                rse_df <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                               rowData=annot,
                                               colData=colData(rse_tx)) %>%
                    as("RangedSummarizedExperiment")
            } else if (region == 'gyrus'){
                rnames <- sapply(strsplit(rownames(colData(rse_tx)), "_"), "[", 1)
                rownames(colData(rse_tx)) <- rnames
                counts <- data.table::fread(tx_file) %>%
                    column_to_rownames("transcript_id")
                colnames(counts) <- sapply(strsplit(colnames(counts), "_"), "[", 1)
                counts <- counts[, rnames] %>% as.matrix
                annot <- data.table::fread(tx_file) %>%
                    column_to_rownames("transcript_id") %>%
                    select(-starts_with("R"))
                rse_df <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                               rowData=annot,
                                               colData=colData(rse_tx)) %>%
                    as("RangedSummarizedExperiment")
            } else {
                counts <- data.table::fread(tx_file) %>%
                    column_to_rownames("transcript_id") %>%
                    select(colnames(rse_tx)) %>% as.matrix
                annot <- data.table::fread(tx_file) %>%
                    column_to_rownames("transcript_id") %>%
                    select(-starts_with("R"))
                rse_df <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                               rowData=annot,
                                               colData=colData(rse_tx)) %>%
                    as("RangedSummarizedExperiment")
                geneid = rownames(rowData(rse_df))
            }
        } else if (exists("rse_exon")){
            rse_df = rse_exon
        } else {
            rse_df = rse_jxn
        }
        keepIndex = which((rse_df$Dx %in% c("Control")) & rse_df$Age > 17 &
                          rse_df$Race %in% c("AA", "CAUC"))
        rse_df = rse_df[, keepIndex]
        rse_df$Sex <- factor(rse_df$Sex)
        rse_df <- merge_rse_metrics(rse_df)
        colData(rse_df)$RIN = sapply(colData(rse_df)$RIN,"[",1)
        rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"), "[", 1)
        pheno = colData(rse_df) %>% as.data.frame %>%
            inner_join(data.table::fread(ancestry),
                       by=c("BrNum"="id", "Race"="group"))
        if(feature == "junctions"){
            datalist1[[region]] <- assays(rse_df)$counts[, pheno$RNum] %>%
                                                as.data.frame %>%
                                                rownames_to_column("V1") %>%
                                                filter(V1 %in% memJuncs()) %>%
                                                column_to_rownames("V1") %>% t %>%
                                                as.data.frame
        }else{
            datalist1[[region]] <- assays(rse_df)$counts[, pheno$RNum] %>% t %>%
                as.data.frame
        }
        datalist2[[region]] <- pheno[, c("Eur", 'Sex', 'Age', 'mitoRate',
                                         'rRNA_rate', 'totalAssignedGene',
                                         'overallMapRate', 'Region', 'Race',
                                         'BrNum', "RNum")]
    }
    samples <- bind_rows(datalist2) %>% as.data.frame %>%
        inner_join(memPROP(), by="RNum") %>%
        mutate_if(is.numeric, scales::rescale)
    if(length(check_dup(samples)) != 0){
        samples <- samples %>% select(-check_dup(samples))
    }
    subjects = samples$BrNum
    if(feature == "junctions"){
        genes = rowData(rse_df)%>% as.data.frame %>%rownames_to_column("V1") %>%
            filter(V1 %in% memJuncs()) %>% column_to_rownames("V1")
    } else {
        genes = rowData(rse_df) %>% as.data.frame
    }
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
    design0 <- model.matrix(~Race*Region, data=x$samples)
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
             'overallMapRate', 'Astro', "Endo", "Macro", "Micro", "Mural",
             "Oligo", "OPC", "Excit", "Inhib", 'Region', 'BrNum')
    adjust_covs = covs[2:15]
    random_effect = covs[17]
    formula <- glue("~ Eur*(1|Region) + ",
                    glue_collapse(adjust_covs, sep = " + ")) %>%
        glue(" + (1|Sex) + (1|{random_effect})")
    return(formula)
}

get_null <- function(){
    covs = c('Sex', 'Age', 'mitoRate', 'rRNA_rate', 'totalAssignedGene',
             'overallMapRate', 'Astro', "Endo", "Macro", "Micro", "Mural",
             "Oligo", "OPC", "Excit", "Inhib", 'Region', 'BrNum')
    adjust_covs = covs[2:15]
    random_effect = covs[17]
    null_form <- glue("~ ", glue_collapse(adjust_covs, sep = " + ")) %>%
        glue(" + (1|Sex) + (1|{random_effect})")
    return(null_form)
}

get_voom <- function(feature, param){
    x <- memPREP(feature)
    vobjDream <- voomWithDreamWeights(x, get_model(), x$samples, BPPARAM=param)
    return(vobjDream)
}
memVOOM <- memoise::memoise(get_voom)

cal_res <- function(feature, param){
### Calculate residuals
    vobj <- memVOOM(feature, param)
    residList <- fitVarPartModel(vobj, get_null(), vobj$targets,
                                 fxn=residuals, BPPARAM=param)
    do.call(rbind,residList)%>% as.data.frame%>%rownames_to_column("Geneid") %>%
        data.table::fwrite(file=paste0(feature, '/residualized_expression.tsv'),
                           sep='\t')
}
memRES <- memoise::memoise(cal_res)

for(feature in c("genes", "transcripts", "exons", "junctions")){
    flush.console()
    dir.create(feature)
    print("Set parallel!")
    param = set_parallel(feature)
    print("Generated residualized expression!")
    memRES(feature, param)
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
