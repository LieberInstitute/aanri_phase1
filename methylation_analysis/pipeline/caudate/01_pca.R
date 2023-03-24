library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)

setwd("out")

# load sd of raw DNAm
i <- 1
load(paste0("chr",i,".rda"))
v <- data.frame(chr=i,start=start(BSobj),sd=sds)
for(i in 2:22){
	cat(i,"\n")
	load(paste0("chr",i,".rda"))
	tmp <- data.frame(chr=i,start=start(BSobj),sd=sds)
	v <- rbind(v,tmp)
}

# get top 1M variable CpG
v <- v[order(v$sd * -1),]
v <- v[1:10^6,]

# DNAm levels for top 1M
i <- 1
tmp <- v[v$chr==i,]
load(paste0("chr",i,".rda"))
BS <- BSobj[is.element(start(BSobj),tmp$start),]
meth <- as.matrix(getMeth(BS))
for(i in 2:22){
	cat(i,"\n")
	load(paste0("chr",i,".rda"))
	tmp <- v[v$chr==i,]
	BS <- BSobj[is.element(start(BSobj),tmp$start),]
	meth <- rbind(meth,as.matrix(getMeth(BS)))
}
id <- colData(BS)$brnum
colnames(meth) <- id

# get ancestry
f_ances <- "/dcs04/lieber/statsgen/shizhong/AANRI/structure/structure_CEU_AFR/structure.out_ancestry_proportion_raceDemo_compare"
ances <- read.table(f_ances,header=T)
idx <- match(id,ances$id)

# residual after regressing out Afr proportion
res <- meth
for(i in 1:nrow(meth)){
	if(! i %% 10000){
		cat(i,"\n")
	}
	model = lm(meth[i,] ~ ances$Afr[idx])
	res[i,] = resid(model)
}

# pca
pc = prcomp(t(res), scale. = T, center = T)
cumsum(summary(pc)$importance[2,1:10])
# PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10
#0.11994 0.16300 0.19064 0.21188 0.22916 0.24498 0.25889 0.27264 0.28560 0.29838

pdf("pca.pdf")
plot(pc)
plot(pc$x[,1], pc$x[,2])
dev.off()
write.csv(pc$x,"pc.csv")

# correlation of pc with ancestry
res <- matrix(NA,nrow=ncol(pc$x),ncol=2)
for(i in 1:ncol(pc$x)){
	tmp <- cor.test(ances$Afr[idx],pc$x[,i])
	res[i,1] <- tmp$estimate
	res[i,2] <- tmp$p.value
}
write.csv(res,"pc_ances_cor.csv")
