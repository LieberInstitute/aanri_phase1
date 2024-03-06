setwd("out_local")
fs <- list.files("./", "_AFR")
res <- c()
for(i in 1:length(fs)){
	if(! i %% 100){
		cat(i,"\n")
	}
	dat <- read.table(fs[i],header=T)
	res <- cbind(res,dat[,2])
}
colnames(res) <- gsub("_AFR","",fs)

id <- read.table("/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC",header=T)
rownames(res) <- id$BrNum[match(dat[,1],id$ID)]
write.csv(res,"../out/local_ancesty.csv")