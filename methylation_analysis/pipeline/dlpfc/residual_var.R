library(data.table)
args = commandArgs(trailingOnly=TRUE)

f_p <- args[1]
f_p_names <- args[2]
chr  <- args[3]
f_out <- args[4]

#f_p <- "/dcs04/lieber/statsgen/shizhong/wgbs/finemap/dlpfc/pheno/chr1/1/p7"
#f_p_names <- "/dcs04/lieber/statsgen/shizhong/wgbs/finemap/dlpfc/pheno/chr1/1/cgnames_p7"
#f_out <- "test";

f_ances <- "/dcs04/lieber/statsgen/shizhong/AANRI/structure/structure_CEU_AFR/structure.out_ancestry_proportion_raceDemo_compare"
f_demo <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC"

# read data
p <- fread(f_p,header=F,data.table=F)
id <- p[,1]
p <- p[,-c(1,2)]
p_names <- read.table(f_p_names,header=F)
ances <- read.table(f_ances,header=T)
demo  <- read.table(f_demo,header=T)

# remove CT snps at CpG sites
f_snp <- paste0("/dcs04/lieber/statsgen/shizhong/AANRI/DEM2/snps_CT/",chr)
snp <- fread(f_snp,header=F,data.table=F)
idx <- is.element(p_names[,1], snp[,1])
p <- p[,!idx]
p_names <- p_names[!idx,1]

# remove low coverage sites
#####################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#change brain region names!!!!!!!!!!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#####################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
load(paste0("/dcs04/lieber/statsgen/shizhong/AANRI/DEM2/dlpfc/aa/coverage5x/",chr,".rda"))
idx <- is.element(p_names,posCovLow)
p <- p[,!idx]
p_names <- p_names[!idx]

# keep AA only
id2 <- intersect(intersect(ances$id[ances$group == "AA"],id),demo$BrNum[demo$Dx == "Control" & demo$Age >= 17])
p <- p[match(id2,id),]

# load pc
pc <- read.csv("./out/pc.csv")

# align samples
ind <- intersect(id2,pc$X)
pc <- pc[match(ind,pc$X),-1]
p <- p[match(ind,id2),]

# get pc to regree out
#pc_ances <- read.csv("./out/pc_ances_cor.csv")
#idx <- which(pc_ances[,3] > 0.05)
#idx <- idx[idx <= 30]

# regress out top5 PC
res <- p
for(i in 1:ncol(p)){
	if(! i %% 100){
		cat(i,"\n")
	}
	d <- as.data.frame(cbind(y=p[,i],pc[,1:5]))
	model = lm(y ~ .,data=d)
	res[,i] = resid(model)
}

# get variance for residual
v <- apply(res,2,sd)

# output
out <- data.frame(chr=chr,pos=p_names,sd=v)
write.table(out,f_out,col.names=F,row.names=F,sep="\t",quote=F)


