suppressMessages({
    library(here)
    library(dplyr)
    library(argparse)
    library(bigstatsr)
})

#### MAIN analysis
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue to be analyzed [default: %default]")
parser$add_argument("-s", "--sge_id", type="integer",
                    help="Chunk to run analysis for")
parser$add_argument("-k", "--k_fold", type="integer", default=5,
                    help="Number of k-folds [default: %default]")
args   <- parser$parse_args()
                                        # Run analysis
feature <- args$feature
tissue  <- args$tissue
sge_id  <- args$sge_id
k_fold  <- args$k_fold

                                        # read covariates
regions <- list("caudate"="caudate", "dentateGyrus"="dg",
                "dlpfc"="dlpfc", "hippocampus"="hippo")
fn <- here("differential_analysis/", tissue,
           "_m/", feature,"/voomSVA.RData")
load(fn)
cov <- v$design %>% as.data.frame %>%
    select(-Intercept, -EA)
sfile   <- here("input/phenotypes/_m",
                paste0(regions[[tissue]],"_phenotypes.csv"))
sp  <- read.csv(sfile, row.names=1) %>%
    select("BrNum", "RNum")
cov <- merge(sp, cov, by=0) %>%
    tibble::column_to_rownames("BrNum") %>%
    select(-"Row.names", -"RNum")

                                        # read gene expression data
outdir  <- paste(tissue, feature, "out", sep="/")
pfile   <- paste0(outdir, "/p", sge_id)
p       <- read.table(pfile, comment.char="", header=TRUE)

                                        # gene meta 
chr   <- gsub("chr","",p[,1])
wind  <- 500000
start <- ifelse(p[,2] - wind < 0, 0, p[,2] - wind)
end   <- p[,3] + wind
genes <- p[,4]

                                        # align p with cov
p   <- t(p[,-(1:4)])
idx <- intersect(rownames(cov),rownames(p))
p   <- p[match(idx,rownames(p)),]
cov <- cov[match(idx,rownames(cov)),]

                                        # keep samples with genotypes
pfam <- paste(tissue,feature,"gwas/chr1_brnum.psam",sep="/")
ind  <- read.table(pfam,comment.char="", header=TRUE)
id   <- intersect(rownames(p),ind[,1])
p    <- p[is.element(rownames(p),id),]
cov  <- cov[is.element(rownames(cov),id),]

                                        # regress out covariates
res <- p
for(i in 1:length(genes)){
    d       <- as.data.frame(cbind(y=p[,i],cov))
	model   <- lm(y ~ .,data=d)
	res[,i] <- resid(model)
}

                                        # fit models
set.seed(20221104)
                                        # randomize sample
N   <- nrow(p)
idx <- sample(1:N)
p   <- p[idx,]
cov <- cov[idx,]
res <- res[idx,]

                                        # k-fold CV
plinkfile0 <- paste(tissue, feature, "gwas/chr", sep="/")
nfold      <- k_fold
folds      <- cut(seq(1,N),breaks=nfold,labels=FALSE)
cor_mat    <- matrix(0,nrow=length(genes),ncol=5)

res_pred   <- res
for(j in 1:length(genes)){
    cat(j,"\n")
                                        # get genotypes
    chr1 <- chr[j]
    p1   <- start[j]
    p2   <- end[j]
    plinkfile <- paste0(plinkfile0,chr1,"_brnum")
    outfile   <- paste0(pfile,"_temp")
    command   <- paste("plink2 --pfile", plinkfile,
                       "--chr", chr1, "--from-bp", p1,
                       "--to-bp", p2, "--export A",
                       "--out", outfile, sep=" ")
    system(command)
                                        # next gene if no snps
    if(! file.exists(paste0(outfile,".raw"))){
        res_pred[,j]  <- NA
        cor_mat[j,]   <- NA
        next		
    }
    geno <- read.table(paste0(outfile,".raw"), header=TRUE)
                                        # Drop SNPs with missing values
    geno <- geno[, apply(geno, 2, function(x) !any(is.na(x)))]

                                        # next gene if no snps
                                        # after quality control
    if(dim(geno)[2] == 6){
        res_pred[,j]  <- NA
        cor_mat[j,]   <- NA
        next
    }
    geno <- as_FBM(as.data.frame(geno[match(rownames(res),geno$FID),-(1:6)]))
    
                                        # next gene if only one snp,
                                        # as big_spLinReg needs >= 2 snps
    if(dim(geno)[2] == 1){
        res_pred[,j]  <- NA
        cor_mat[j,]   <- NA
        next		
    }
                                        # model and prediction
    betas1 <- list()
    for(i in 1:nfold){
        ind.test  <- which(folds==i,arr.ind=TRUE)
        ind.train <- setdiff(rows_along(geno), ind.test)
                                        # model based on expression residuals
        mod  <- big_spLinReg(geno, res[ind.train,j],
                             ind.train = ind.train,
                             alphas = seq(0.05,1,0.05), K = 4)
        pred <- predict(mod, geno, ind.test)
        betas1[[i]] <- data.frame(SNP=attr(mod, "ind.col"),
                                  beta = summary(mod, best.only=T)$beta[[1]])
        cor_mat[j,i] <- cor(pred,res[ind.test,j])
    }
                                        # residuals
    weights_df    <- apply(purrr::reduce(betas1, full_join, by="SNP") %>%
                           tibble::column_to_rownames("SNP"), 1,
                           function(x) mean(x, na.rm=TRUE))
    X  <- as_FBM(as.data.frame(geno[, as.numeric(names(weights_df))]))
    res_pred[,j]  <- big_prodVec(X, weights_df)
}

                                        # save output files
write.csv(cbind(genes, cor_mat),
          paste0(pfile,"_cor_res.csv"))
colnames(res_pred)  <- genes
write.csv(res_pred, paste0(pfile,"_res_pred.csv"))
                                        # clean up
outfile <- paste0(pfile,"_temp.raw")
system(paste("rm", outfile, sep=" "))
outfile <- paste0(pfile,"_temp.log")
system(paste("rm", outfile, sep=" "))

#### Reproducibility information
if(sge_id == 1){
    Sys.time()
    proc.time()
    options(width = 120)
    sessioninfo::session_info()
}
