library(bigstatsr)
args = commandArgs(trailingOnly=TRUE)

# read gene expression data and covariates
pfile <- args[1]
covfile <- "/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/eqtl_analysis/caudate/genes/covariates/_m/genes.combined_covariates.txt"
p <- read.table(pfile,comment.char="", header=T)
cov <- read.table(covfile,comment.char="", header=T)

# gene meta 
chr <- gsub("chr","",p[,1])
wind <- 500000
start <- ifelse(p[,2] - wind < 0, 0, p[,2] - wind)
end <- p[,3] + wind
genes <- p[,4]

# align p with cov
p <- t(p[,-(1:4)])
rownames(cov) <- cov[,1]
cov <- t(cov[,-1])
idx <- intersect(rownames(cov),rownames(p))
p <- p[match(idx,rownames(p)),]
cov <- cov[match(idx,rownames(cov)),]

# keep samples with genotypes
pfam <- "./gwas/chr1_brnum.psam"
ind <- read.table(pfam,comment.char="", header=T)
id <- intersect(rownames(p),ind[,1])
p <- p[is.element(rownames(p),id),]
cov <- cov[is.element(rownames(cov),id),]

# regress out covariates
res <- p
for(i in 1:length(genes)){
	d <- as.data.frame(cbind(y=p[,i],cov))
	model = lm(y ~ .,data=d)
	res[,i] = resid(model)
}

# fit models
set.seed(2022)

plinkfile0 <- "./gwas/chr"
res_pred <- res
for(j in 1:length(genes)){
	cat(j,"\n")
	# get genotypes
	chr1 <- chr[j]
	p1 <- start[j]
	p2 <- end[j]
	plinkfile <- paste0(plinkfile0,chr1,"_brnum")
	outfile <- paste0(pfile,"_temp")
	command <- paste("plink2 --pfile", plinkfile,"--chr",chr1,"--from-bp",p1,"--to-bp",p2,"--export A", "--out", outfile,sep=" ")
	system(command)
	# next gene if no snps
	if(! file.exists(paste0(outfile,".raw"))){
		res_pred[,j] <- NA
		next		
	}
	geno <- read.table(paste0(outfile,".raw"),header=T)
	geno <-as_FBM(as.data.frame(geno[match(rownames(res),geno$IID),-(1:6)]))
	# next gene if only one snp, as big_spLinReg needs >= 2 snps
	if(dim(geno)[2] == 1){
		res_pred[,j] <- NA
		next		
	}
	# model based on expression residuals
	mod <- big_spLinReg(geno, res[,j], K = 4,alphas = seq(0.1,1,0.1))
	res_pred[,j] <- predict(mod, geno)
}

# save output files
outfile <- paste0(pfile,"_res_pred.csv")
write.csv(cbind(genes,res_pred),outfile)

# clean up
outfile <- paste0(pfile,"_temp.raw")
system(paste("rm", outfile,sep=" "))
outfile <- paste0(pfile,"_temp.log")
system(paste("rm", outfile,sep=" "))
