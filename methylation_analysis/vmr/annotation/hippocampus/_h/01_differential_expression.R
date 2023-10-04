## Differential expression with limma-voom pipeline
suppressMessages({
    library(SummarizedExperiment)
    library(argparse)
    library(dplyr)
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

get_vmr_resid <- function(){
    vmr_file <- "../../../_m/aaOnly/hippocampus_vmr_residual.txt"
    return(data.table::fread(vmr_file))
}
memRES <- memoise::memoise(get_vmr_resid)

get_de <- function(feature){
    feat_lt <- list("genes"="Gene", "transcripts"="Transcript",
                    "exons"="Exon", "junctions"="Junction")
    de_file <- here::here("differential_analysis/tissue_comparison",
                          "summary_table/_m",
                          "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    return(data.table::fread(de_file) |>
           filter(Type == feat_lt[[feature]], Tissue == "Hippocampus"))
}
memDE <- memoise::memoise(get_de)

get_sig_de <- function(feature){
    return(memDE(feature) |> filter(lfsr < 0.05))
}
memSIG <- memoise::memoise(get_sig_de)

annot_vmrs <- function(symbol){
    vmr_file <- "../../_m/hippocampus/vmr_annotation.tsv"
    return(data.table::fread(vmr_file) |>
           select(vmr_id, seqnames, start, end, annot.symbol) |>
           filter(annot.symbol == symbol) |> distinct())
}
memVMRS <- memoise::memoise(annot_vmrs)

select_vmrs <- function(feature, geneid){
    symbol <- filter(memSIG(feature), Effect == geneid)$Symbol
    if(symbol == ""){
        return(NULL)
    }
    if(dim(memVMRS(symbol))[1] == 0){
        return(NULL)
    }
    x      <- memPREP(feature)
    vmr_id <- memVMRS(symbol)$vmr_id
    vmr_df <- memRES() |> select("V1", any_of(vmr_id)) |>
        rename("BrNum"="V1")
    return(x$samples |> as.data.frame() |> select(RNum, BrNum) |>
           inner_join(vmr_df, by="BrNum") |> select(-"BrNum"))
}

prep_data <- function(feature){
    ancestry = here::here("input/ancestry_structure",
                          "structure.out_ancestry_proportion_raceDemo_compare")
    counts_lt = list("genes"=paste0(here::here(), "/input/counts/_m/",
                                    "hippo_brainseq_phase2_hg38_",
                                    "rseGene_merged_n447.rda"),
                     "transcripts"=paste0(here::here(), "/input/counts/_m/",
                                    "hippo_brainseq_phase2_hg38_",
                                    "rseTx_merged_n447.rda"),
                     "exons"=paste0(here::here(), "/input/counts/_m/",
                                    "hippo_brainseq_phase2_hg38_",
                                    "rseExon_merged_n447.rda"),
                     "junctions"=paste0(here::here(), "/input/counts/_m/",
                                    "hippo_brainseq_phase2_hg38_",
                                    "rseJxn_merged_n447.rda"))
    tx_file = here::here("input/counts/_m","transcripts_counts",
                         "hippo_counts.txt")
    load(counts_lt[[feature]], verbose=TRUE)
    if(exists("rse_gene")){
        rse_df  <- rse_gene
    } else if (exists("rse_tx")){
        rnames <- sapply(strsplit(rownames(colData(rse_tx)), "_"), "[", 1)
        counts <- data.table::fread(tx_file) |>
            tibble::column_to_rownames("transcript_id") 
        colnames(counts) <- sapply(strsplit(colnames(counts), "_"), "[", 1)
        counts <- counts |> select(any_of(rnames)) |> as.matrix()
        annot <- rowData(rse_tx) |> as.data.frame() |> 
            filter(transcript_id %in% rownames(counts)) |> as.matrix()
        keepIdx <- which(rnames %in% colnames(counts))
        rse_tx  <- rse_tx[, keepIdx]
        rse_df <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                       rowData=annot, colData=colData(rse_tx)) |>
            as("RangedSummarizedExperiment")
        colnames(rse_df) <- rnames[keepIdx]
    } else if (exists("rse_exon")){
        rse_df  <- rse_exon
    } else {
        rse_df  <- rse_jxn
    }
    keepIndex   <- which((rse_df$Dx %in% c("Control")) & 
                         rse_df$Age > 17 & rse_df$Race %in% c("AA"))
    rse_df      <- rse_df[, keepIndex]
    rse_df$Dx   <- factor(rse_df$Dx, levels = c("Control", "Schizo"))
    rse_df$Sex  <- factor(rse_df$Sex)
    rse_df      <- merge_rse_metrics(rse_df)
    colData(rse_df)$RIN       <- sapply(colData(rse_df)$RIN,"[",1)
    rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"),
                                        "[", 1)
                                        # Generate phenotypes
    pheno = colData(rse_df) |> as.data.frame() |> 
        inner_join(data.table::fread(ancestry),
                   by=c("BrNum"="id", "Race"="group")) |>
        filter(BrNum %in% memRES()$V1)
                                        # Generate DGE list
    x <- edgeR::DGEList(counts=assays(rse_df)$counts[, pheno$RNum], 
                        genes=rowData(rse_df), samples=pheno)
    x <- x[memDE(feature)$Effect, , keep.lib.sizes=FALSE]
    x <- edgeR::calcNormFactors(x, method="TMM")
    return(x)
}
memPREP <- memoise::memoise(prep_data)

gen_model <- function(feature){
    x <- memPREP(feature)
    # Design matrix
    mod = model.matrix(~Sex + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + overallMapRate, data = x$samples)
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Load qSV
    qsv_file = here::here("input/phenotypes/_m/qSV_hippo.csv")
    modQsva <- mod |> as.data.frame() |> tibble::rownames_to_column("RNum") |>
        inner_join(data.table::fread(qsv_file), by=c("RNum"="V1")) |> 
        rename_all(list(~stringr::str_replace_all(., 'PC', 'qPC'))) |>
        tibble::column_to_rownames("RNum") |> as.matrix()
    return(modQsva)
}
memGEN <- memoise::memoise(gen_model)

qSV_model <- function(feature){
    x <- memPREP(feature)
    # Design matrix
    mod = model.matrix(~Eur + Sex + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + overallMapRate, data = x$samples)
    colnames(mod) <- gsub("Eur", "EA", colnames(mod))
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Load qSV
    qsv_file = here::here("input/phenotypes/_m/qSV_hippo.csv")
    modQsva <- mod |> as.data.frame() |> tibble::rownames_to_column("RNum") |>
        inner_join(data.table::fread(qsv_file), by=c("RNum"="V1")) |> 
        rename_all(list(~stringr::str_replace_all(., 'PC', 'qPC'))) |>
        tibble::column_to_rownames("RNum") |> as.matrix()
    return(modQsva)
}
memQSV <- memoise::memoise(qSV_model)

vmr_model <- function(feature, geneid){
    x <- memPREP(feature)
    # Design matrix
    mod = model.matrix(~Eur + Sex + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + overallMapRate, data = x$samples)
    colnames(mod) <- gsub("Eur", "EA", colnames(mod))
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Load qSV
    qsv_file <- here::here("input/phenotypes/_m/qSV_hippo.csv")
    modQsva  <- mod |> as.data.frame() |>
        tibble::rownames_to_column("RNum") |>
        inner_join(data.table::fread(qsv_file), by=c("RNum"="V1")) |> 
        rename_all(list(~stringr::str_replace_all(., 'PC', 'qPC'))) |>
        inner_join(select_vmrs(feature, geneid), by="RNum") |>
        tibble::column_to_rownames("RNum") |> as.matrix()
    return(modQsva)
}
memVMR <- memoise::memoise(vmr_model)

vmr_model_gen <- function(feature, geneid){
    x <- memPREP(feature)
    # Design matrix
    mod = model.matrix(~Sex + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + overallMapRate, data = x$samples)
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Load qSV
    qsv_file <- here::here("input/phenotypes/_m/qSV_hippo.csv")
    modQsva  <- mod |> as.data.frame() |>
        tibble::rownames_to_column("RNum") |>
        inner_join(data.table::fread(qsv_file), by=c("RNum"="V1")) |> 
        rename_all(list(~stringr::str_replace_all(., 'PC', 'qPC'))) |>
        inner_join(select_vmrs(feature, geneid), by="RNum") |>
        tibble::column_to_rownames("RNum") |> as.matrix()
    return(modQsva)
}
memVMR2 <- memoise::memoise(vmr_model_gen)

get_norm <- function(feature, geneid, fnc, VMR){
    ### Preform voom
    x <- memPREP(feature)
    if(VMR){
        modQsva <- fnc(feature, geneid)
    } else {
        modQsva <- fnc(feature)
    }
    y <- edgeR::cpm(x[geneid, rownames(modQsva)], log=TRUE)
    return(y)
}
memNORM <- memoise::memoise(get_norm)

fit_model <- function(feature, geneid, fnc, VMR){
    y       <- memNORM(feature, geneid, fnc, VMR)
    if(VMR){
        modQsva <- fnc(feature, geneid)
    } else {
        modQsva <- fnc(feature)
    }
    fit <- limma::lmFit(y, modQsva)
    return(fit)
}
memFIT <- memoise::memoise(fit_model)

get_partial_r2 <- function(feature, feature_id, GEN){
                                        # Load normalized expression
    y    <- memNORM(feature, feature_id, memQSV, FALSE) |>
        t() |> as.data.frame() |>
        tibble::rownames_to_column("RNum")
    if(GEN){
        label = "general"
    } else {
        label = "ancestry"
    }
                                        # Load covariates
    if(GEN){
        ## reduced model
        mod1 <- memGEN(feature) |> as.data.frame() |>
            select(-Intercept) |>
            tibble::rownames_to_column("RNum")
        ## full model
        mod2 <- memQSV(feature) |> as.data.frame() |>
            select(-Intercept) |>
            tibble::rownames_to_column("RNum")
    } else {
        ## reduced model
        mod1 <- memVMR2(feature, feature_id) |>
            as.data.frame() |> select(-Intercept) |>
            tibble::rownames_to_column("RNum")
        ## full model
        mod2 <- memVMR(feature, feature_id) |>
            as.data.frame() |> select(-Intercept) |>
            tibble::rownames_to_column("RNum")
    }
                                        # Combine data
    df1  <- inner_join(y, mod1, by="RNum") |>
        tibble::column_to_rownames("RNum")
    df2  <- inner_join(y, mod2, by="RNum") |>
        tibble::column_to_rownames("RNum")
                                        # Generate linear model
    model1 <- paste0(feature_id, " ~ ",
                     paste(colnames(df1)[-1], collapse=" + "))
    model2 <- paste0(feature_id, " ~ ",
                     paste(colnames(df2)[-1], collapse=" + "))
                                        # Partial R2
    reduced <- anova(lm(model1, data=df1))
    full    <- anova(lm(model2, data=df2))
    p1      <- reduced["Residuals", "Sum Sq"]
    p2      <- full["Residuals", "Sum Sq"]
    partial_r2 <- (p1 - p2) / p1
    return(data.frame("feature_id"=feature_id,
                      "partial_r2"=partial_r2,
                      "full_r2"=p2, "reduced_r2"=p1,
                      "Model"=label))
}

## extract_de_loop <- function(feature, feature_id, VMR){
##     efit   <- memFIT(feature, feature_id, VMR)
##     bhat   <- efit$coefficients[, 2]
##     shat   <- efit$stdev.unscaled[, 2] / sqrt(dim(efit$design)[1])
##     return(data.frame("feature_id"=feature_id, "bhat"=bhat, "shat"=shat[1]))
## }

#### MAIN analysis
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-c", "--chunks", type="integer", default=20,
                    help="Chunks to analysis data [default: %default]")
parser$add_argument("-t", "--sge_id", type="integer",
                    help="SGE ID to select geneid")
args   <- parser$parse_args()
                                        # Prep inputs
feature <- args$feature; n <- args$chunks; ID <- args$sge_id
if(!dir.exists(feature)){ dir.create(feature) }

                                        # Chunk data for loop
ids <- memSIG(feature)$Effect |> sort()
ps  <- cut(seq(1,length(ids)), breaks=n, labels=FALSE)
idx <- which(ps == ID, arr.ind=TRUE)
select_ids <- ids[idx]

                                        # DE by gene
for(feature_id in select_ids){
    if(!is.null(select_vmrs(feature, feature_id))){
                                        ## # no VMR
        ## outfile1 <- paste0(feature, "/diffExpr_model.",
        ##                    feature_id, ".main")
        ## extract_de_loop(feature, feature_id, FALSE) |>
        ##     data.table::fwrite(outfile1, sep='\t')
        ##                                 # VMR
        ## outfile2 <- paste0(feature, "/diffExpr_model.",
        ##                    feature_id, ".vmr")
        ## extract_de_loop(feature, feature_id, TRUE) |>
        ##     data.table::fwrite(outfile2, sep='\t')
                                        # Calculate partial R2
        outfile1 <- paste0(feature, "/partial_r2.",
                           feature_id, ".general.r2")
        get_partial_r2(feature, feature_id, TRUE) |>
            data.table::fwrite(outfile1, sep='\t')
        outfile2 <- paste0(feature, "/partial_r2.",
                           feature_id, ".ancestry.r2")
        get_partial_r2(feature, feature_id, FALSE) |>
            data.table::fwrite(outfile2, sep='\t')
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
