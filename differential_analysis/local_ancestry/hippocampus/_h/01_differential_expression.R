## Differential expression with limma-voom pipeline
suppressMessages({
    library(SummarizedExperiment)
    library(argparse)
    library(dplyr)
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

prep_data <- function(feature){
    counts_lt = list("genes"=here("input/counts/_m",
                                  paste0("hippo_brainseq_phase2_",
                                         "hg38_rseGene_merged_n447.rda")),
                     "transcripts"=here("input/counts/_m",
                                        paste0("hippo_brainseq_phase2_",
                                               "hg38_rseTx_merged_n447.rda")),
                     "exons"=here("input/counts/_m",
                                  paste0("hippo_brainseq_phase2_",
                                         "hg38_rseExon_merged_n447.rda")),
                     "junctions"=here("input/counts/_m",
                                      paste0("hippo_brainseq_phase2_",
                                             "hg38_rseJxn_merged_n447.rda")))
    tx_file = here("input/counts/_m/transcripts_counts/hippo_counts.txt")
    load(counts_lt[[feature]], verbose=TRUE)
    if(exists("rse_gene")){
        rse_df = rse_gene
    } else if (exists("rse_tx")){
        counts <- data.table::fread(tx_file) |>
            tibble::column_to_rownames("transcript_id") |>
            select(any_of(colnames(rse_tx))) |> as.matrix()
        annot <- rowData(rse_tx) |> as.data.frame() |> 
            filter(transcript_id %in% rownames(counts)) |> as.matrix()
        keepIdx <- which(colnames(rse_tx) %in% colnames(counts))
        rse_tx  <- rse_tx[, keepIdx]
        rse_df <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                       rowData=annot, colData=colData(rse_tx)) |>
            as("RangedSummarizedExperiment")
    } else if (exists("rse_exon")){
        rse_df = rse_exon
    } else {
        rse_df = rse_jxn
    }
    keepIndex = which((rse_df$Dx %in% c("Control")) & 
                      rse_df$Age > 17 & rse_df$Race %in% c("AA"))
    rse_df = rse_df[, keepIndex]
    rse_df$Dx = factor(rse_df$Dx, levels = c("Control", "Schizo"))
    rse_df$Sex <- factor(rse_df$Sex)
    rse_df <- merge_rse_metrics(rse_df)
    colData(rse_df)$RIN = sapply(colData(rse_df)$RIN,"[",1)
    rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"),
                                        "[", 1)
                                        # Generate phenotypes
    pheno = colData(rse_df) |> as.data.frame()
                                        # Generate DGE list
    x <- edgeR::DGEList(counts=assays(rse_df)$counts[, pheno$RNum], 
                        genes=rowData(rse_df), samples=pheno)
    keep.x <- edgeR::filterByExpr(x)
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    x <- edgeR::calcNormFactors(x, method="TMM")
    return(x)
}
memPREP <- memoise::memoise(prep_data)

get_fam <- function(){
    plink_file <- here::here("input/genotypes/_m/TOPMed_LIBD_AA.psam")
    return(data.table::fread(plink_file) |>
           select(-SEX, -PHENO1) |> rename("FID"="#FID"))
}
memFAM <- memoise::memoise(get_fam)

get_local_ancestry <- function(feature, feature_id){
    afile = paste0("../../_m/hippocampus/",feature,"/ancestry/",
                   feature_id,"_AFR")
    return(data.table::fread(afile) |>
           inner_join(memFAM(), by=c("ID"="IID")) |>
           select(FID, AFR))
}
memANCESTRY <- memoise::memoise(get_local_ancestry)

filter_genotypes <- function(feature, feature_id){
    x <- memPREP(feature)
                                            # Subset on genotypes
    new_pheno <- x$samples |>
        inner_join(get_local_ancestry(feature, feature_id), by=c("BrNum"="FID"))
    x <- x[, new_pheno$RNum, keep.lib.size=TRUE]
    x$samples <- new_pheno
    rownames(x$samples) <- new_pheno$RNum
    return(x)
}
memFILTER <- memoise::memoise(filter_genotypes)

eqtl_model <- function(feature, feature_id){
    x <- memFILTER(feature, feature_id)
    # Design matrix
    mod = model.matrix(~AFR + Sex + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + overallMapRate, data = x$samples)
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Load qSV
    qsv_file = here("input/phenotypes/_m/qSV_hippo.csv")
    modQsva <- mod |> as.data.frame() |> tibble::rownames_to_column("RNum") |>
        inner_join(data.table::fread(qsv_file), by=c("RNum"="V1")) |> 
        rename_all(list(~stringr::str_replace_all(., 'PC', 'qPC'))) |>
        tibble::column_to_rownames("RNum") |> as.matrix()
    return(modQsva)
}
memMODEL <- memoise::memoise(eqtl_model)

get_norm <- function(feature, feature_id){
    ### Preform voom
    x <- memFILTER(feature, feature_id)
    modQsva <- memMODEL(feature, feature_id)
    y <- edgeR::cpm(x[feature_id, rownames(modQsva)], log=TRUE)
    return(y)
}
memNORM <- memoise::memoise(get_norm)

fit_model <- function(feature, feature_id){
    y       <- memNORM(feature, feature_id)
    modQsva <- memMODEL(feature, feature_id)
    fit <- limma::lmFit(y, modQsva)
    return(fit)
}
memFIT <- memoise::memoise(fit_model)

extract_de_loop <- function(feature, ID){
    efit   <- memFIT(feature, feature_id)
    bhat   <- efit$coefficients[, 2] ## local ancestry
    shat   <- efit$stdev.unscaled[, 2] / sqrt(dim(efit$design)[1])
    return(data.frame("feature_id"=feature_id, "bhat"=bhat, "shat"=shat[1]))
}

#### MAIN analysis
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-c", "--chunks", type="integer", default=200,
                    help="Chunks to analysis data [default: %default]")
parser$add_argument("-t", "--sge_id", type="integer",
                    help="SGE ID to select feature_id")
args <- parser$parse_args()
                                        # Prep inputs
feature <- args$feature; n <- args$chunks; ID <- args$sge_id
                                        # Chunk data for loop
ids <- gsub("_AFR", "",
            list.files(paste0("../../_m/hippocampus/",feature,"/ancestry"),
                       pattern="\\_AFR$")) |> sort()
ps  <- cut(seq(1,length(ids)),breaks=n,labels=FALSE)
idx <- which(ps == ID,arr.ind=TRUE)
select_ids <- ids[idx]
                                        # By gene, local ancestry
for(feature_id in select_ids){
    outfile <- paste0(feature, "/diffExpr_model.",feature_id,".ancestry")
    extract_de_loop(feature, feature_id) |>
        data.table::fwrite(outfile, sep='\t')
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
