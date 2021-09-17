## Calculate the partial R2 using RaFFE SNPs and individually

suppressPackageStartupMessages({
    library("argparse")
    library("tidyverse")
})

get_ml_summary <- function(fn){
    ml_df = data.table::fread(fn) %>% mutate_at("fold", as.character) %>%
        select(fold, n_features, n_redundant, starts_with("test_score_r2")) %>%
        pivot_longer(-fold) %>% group_by(name) %>%
        summarise(Mean=mean(value), Median=median(value), Std=sd(value),
                  .groups = "keep") %>%
        as.data.frame %>% column_to_rownames("name")
    return(ml_df)
}

get_snps_ml <- function(fn1, fn2){
    snps = data.table::fread(fn2) %>%
        rename('Geneid'='V1', 'Fold'='V2', 'Rank'='V3') %>%
        group_by(Fold, Geneid) %>% slice(1) %>%
        arrange(Fold, desc(Rank)) %>% as.data.frame %>%
        pivot_wider(names_from=Fold, values_from=Rank) %>%
        mutate(median_all = apply(., 1, median)) %>%
        arrange(median_all) %>% mutate(Rank = rank(median_all)) %>%
        filter(Rank <= get_ml_summary(fn1)["n_features", "Median"]) %>%
        select(Geneid)
    return(snps)
}

get_pheno <- function(tissue, target, qsv_dir, pheno_file){
    qSV_lt = list("Caudate"="qSV_caudate.csv", "DLPFC"="qSV_dlpfc.csv",
                  "Hippocampus"="qSV_hippo.csv", "Dentate Gyrus"="qSV_dg.csv")
    ancestry = data.table::fread(target)
    qSV = data.table::fread(paste(qsv_dir, qSV_lt[[tissue]], sep="/")) %>%
        rename_all(list(~str_replace_all(., "PC", "qPC")))
    pheno = data.table::fread(pheno_file) %>%
        inner_join(ancestry, by=c("BrNum"="id")) %>%
        inner_join(qSV, by="V1")
    return(pheno)
}

memPHENO <- memoise::memoise(get_pheno)

tissue_map <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus", "dentateGyrus"="Dentate Gyrus")[[tissue]])
}

remove_aliased <- function(model, dfx){
    test = alias(lm(model, data=dfx))$Complete
    for(ii in rownames(test)){
        model = str_remove(model, paste0(ii, " +"))
    }
    return(model)
}

calculate_idv_partial <- function(tissue, gname, target, qsv_dir, pheno_file,
                                  snp_dir){
                                        # Get data
    gene_id = gsub("_", ".", gname)
    dir.create(paste(tissue, gene_id, sep='/'))
    fn = paste(snp_dir, tissue, gname, "snps_onehot.csv", sep="/")
    df = memPHENO(tissue_map(tissue), target, qsv_dir, pheno_file) %>%
        inner_join(data.table::fread(fn), by="V1") %>%
        column_to_rownames("V1") %>% rename_with(~gsub(":", "_", .),
                                                 starts_with('chr'))
                                        # Model 1
    model1 = paste(paste0(gene_id, "~ Eur"), "Sex + Age + mitoRate + rRNA_rate",
                   "overallMapRate", paste(colnames(df)[grep("^qPC", colnames(df))],
                                           collapse=" + "), sep=" + ")
    reduced = anova(lm(model1, data=df))
                                        # Loop through SNPs for model 2
    snp_lt <- c(); partial1 <-  c(); partial2 <-  c(); partial_r2 <-  c()
    for(col in colnames(df)[(grep("^chr", colnames(df)))]){
        snp_lt = c(snp_lt, col)
        model2 = paste(model1, col, sep=" + ")
        full = anova(lm(model2, data=df))
        p1 = reduced["Residuals", "Sum Sq"]
        p2 = full["Residuals", "Sum Sq"]
        partial1 = c(partial1, p1); partial2 = c(partial2, p2)
        partial_r2 = c(partial_r2, (p1 - p2) / p1)
    }
    dt = data.frame("SNP"=snp_lt, "Partial_R2"=partial_r2,
                    "Full_R2"=partial2, "Reduced_R2"=partial1)
    dt["Tissue"] = tissue_map(tissue)
    dt["Geneid"] = gene_id
    dt %>% data.table::fwrite(paste(tissue, gene_id, "individual_partial_r2.tsv", sep='/'),
                              sep='\t')
}

calculate_rf_partial <- function(tissue, target, qsv_dir, pheno_file, ml_dir,
                                 snp_dir){
                                        # Get data
    snp_len <- c(); genes <- c(); partial1 <- c(); partial2 <- c(); partial_r2 <- c()
    ## PATH for onehot files
    glist <- list.files(paste0(snp_dir, tissue, '/'))
    for(gname in glist[!grepl("*\\.log", glist)]){
                                        # Clean gene names
        gene_id = gsub("_", ".", gname)
        #print(gene_id)
                                        # Get data
        fn = paste(snp_dir, tissue, gname, "snps_onehot.csv", sep="/")
        if(file.exists(fn)){
            df = memPHENO(tissue_map(tissue), target, qsv_dir, pheno_file) %>%
                inner_join(data.table::fread(fn), by="V1") %>%
                column_to_rownames("V1") %>%
                rename_with(~gsub(":", "_", .), starts_with('chr'))
                                        # Model 1
            model1 = paste(paste0(gene_id,"~ Eur"),"Sex + Age + mitoRate + rRNA_rate + overallMapRate",
                           paste(colnames(df)[grep("^qPC", colnames(df))], collapse=" + "), sep=" + ")
            reduced = anova(lm(model1, data=df))
                                        # Get RF predictive SNPs
            fn1 = paste(ml_dir, tissue, gname, "rf_10folds.txt", sep="/")
            fn2 = paste(ml_dir, tissue, gname, "rank_features.txt", sep="/")
            snps = get_snps_ml(fn1, fn2) %>% mutate("SNPs"=gsub(":", "_", Geneid))
                                        # Model 2
            if(length(snps$SNPs) == 0){
                print("There are no SNPs")
                next
            } else {
                model2 = paste(model1, paste(snps$SNPs, collapse=" + "), sep=" + ")
                full = anova(lm(model2, data=df))
                p1 = reduced["Residuals", "Sum Sq"]
                p2 = full["Residuals", "Sum Sq"]
                genes = c(genes, gene_id)
                snp_len = c(snp_len, length(snps$SNPs))
                partial1 = c(partial1, p1); partial2 = c(partial2, p2)
                partial_r2 = c(partial_r2, (p1 - p2) / p1)
            }
        }
    }
    dt = data.frame("Geneid"=genes, "N_Features"= snp_len, "Partial_R2"=partial_r2,
                    "Full_R2"=partial2, "Reduced_R2"=partial1)
    dt["Tissue"] = tissue
    dt %>% data.table::fwrite(paste(tissue, "raffe_partial_r2.tsv", sep='/'), sep='\t')
}


calculate_enet_partial <- function(tissue, target, qsv_dir, pheno_file, ml_dir,
                                   snp_dir){
                                        # Get data
    snp_len <- c(); genes <- c(); partial1 <- c(); partial2 <- c(); partial_r2 <- c()
    glist <- list.files(paste0(snp_dir, tissue, '/'))
    for(gname in glist[!grepl("*\\.log", glist)]){
                                        # Clean gene names
        gene_id = gsub("_", ".", gname)
        #print(gene_id)
                                        # Get data
        fn = paste(snp_dir, tissue, gname, "snps_onehot.csv", sep="/")
        df = memPHENO(tissue_map(tissue), target, qsv_dir, pheno_file) %>%
            inner_join(data.table::fread(fn), by="V1") %>%
            column_to_rownames("V1") %>%
            rename_with(~gsub(":", "_", .), starts_with('chr'))
                                        # Model 1
        model1 = paste(paste0(gene_id, "~ Eur"), "Sex + Age + mitoRate + rRNA_rate + overallMapRate",
                       paste(colnames(df)[grep("^qPC", colnames(df))], collapse=" + "), sep=" + ")
        reduced = anova(lm(model1, data=df))
                                        # Get dRFE Elastic Net predictive SNPs
        fn1 = paste(ml_dir, tissue, gname, "enet_rfe_10folds.txt", sep="/")
        fn2 = paste(ml_dir, tissue, gname, "rank_features.txt", sep="/")
        snps = get_snps_ml(fn1, fn2) %>% mutate("SNPs"=gsub(":", "_", Geneid))
                                        # Model 2
        if(length(snps$SNPs) == 0){
            print(gene_id)
            print("There are no SNPs")
            next
        } else {
            model2 = paste(model1, paste(snps$SNPs, collapse=" + "), sep=" + ")
            full = anova(lm(model2, data=df))
            p1 = reduced["Residuals", "Sum Sq"]
            p2 = full["Residuals", "Sum Sq"]
            genes = c(genes, gene_id)
            snp_len = c(snp_len, length(snps$SNPs))
            partial1 = c(partial1, p1); partial2 = c(partial2, p2)
            partial_r2 = c(partial_r2, (p1 - p2) / p1)
        }
    }
    dt = data.frame("Geneid"=genes, "N_Features"= snp_len, "Partial_R2"=partial_r2,
                    "Full_R2"=partial2, "Reduced_R2"=partial1)
    dt["Tissue"] = tissue
    dt %>% data.table::fwrite(paste(tissue, "enet_partial_r2.tsv", sep='/'), sep='\t')
}

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-r", "--rf", action="store_true", default=FALSE,
                    help="Running with dRFE random forest SNPs [default: %(default)s]")
parser$add_argument("-e", "--enet", action="store_true", default=FALSE,
                    help="Running with dRFE elastic net SNPs [default: %(default)s]")
parser$add_argument("-g", "--gname", type="character",
                    help="Name of gene to analyze [required]")
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue for brain analysis [default: %default]")
parser$add_argument("-d", "--snp_dir", type="character", default="../../_m/",
                    help="PATH to onehot encoded SNPs [default: %default]")
parser$add_argument("-x", "--target", type="character",
                    help="File with target for partial R2 (Ex: ancestry) [required]")
parser$add_argument("-S", "--qsv_path", type="character",
                    help="PATH to qSVs [required]")
parser$add_argument("-p", "--pheno_file", type="character",
                    help="PATH to phenotypes [required]")
parser$add_argument("-M", "--ml_dir", type="character",
                    help="PATH to dRFE results [required]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

if(length(args$target) == 0){
    print("Missing target!")
} else if (length(args$qsv_path) == 0){
    print("Missing qSV directory!")
} else if (length(args$pheno_file) == 0){
    print("Missing phenotypes file!")
} else {
    if(args$rf){
        if (length(args$ml_dir) == 0){
            print("Missing dRFE directory!")
        } else {
            calculate_rf_partial(args$tissue, args$target, args$qsv_dir,
                                 args$pheno_file, args$ml_dir, args$snp_dir)
        }
    }else if(args$enet){
        if (length(args$ml_dir) == 0){
            print("Missing dRFE directory!")
        } else {
            calculate_enet_partial(args$tissue, args$target, args$qsv_dir,
                                 args$pheno_file, args$ml_dir, args$snp_dir)
        }
    }else{
        if(length(args$gname) == 0){
            print("Missing gname!")
        } else {
            calculate_idv_partial(args$tissue, args$gname)
        }
    }
}
