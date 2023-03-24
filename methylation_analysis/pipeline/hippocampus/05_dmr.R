library(rGREAT)

setwd("out")

# load DNAm levels per region and all vmrs
load("chr1_vmr.rda")
meth2 <- meth
vmrs2 <- vmrs
for(i in 2:22){
	load(paste0("chr",i,"_vmr.rda"))
	meth2 <- rbind(meth2,meth)
	vmrs2 <- rbind(vmrs2,vmrs)
}
p <- t(meth2)
ind <- rownames(p)

# replace missing data by col means
for(i in 1:ncol(p)){
	p[is.na(p[,i]),i] <- mean(p[,i],na.rm=T)
}

# read covariates
f_ances <- "/dcs04/lieber/statsgen/shizhong/AANRI/structure/structure_CEU_AFR/structure.out_ancestry_proportion_raceDemo_compare"
f_demo <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC"
f_pc <- "pc.csv"
ances <- read.table(f_ances,header=T)
demo  <- read.table(f_demo,header=T)
ances <- ances[match(ind,ances$id),]
demo <- demo[match(ind,demo$BrNum),]
pc <- read.csv(f_pc)
pc <- pc[match(ind,pc$X),-1]

# differential DNAm for Afr proportion
res <- matrix(NA,nrow=ncol(p),ncol=4)
colnames(res) <- c("beta","se", "t", "p")
for(i in 1:ncol(p)){
	if(! i %% 100){
		cat(i,"\n")
	}
	d <- as.data.frame(cbind(y=p[,i],ances$Afr,pc[,1:5],demo$Age,as.factor(demo$Sex)))
	model = lm(y ~ .,data=d)
	res[i,] = summary(model)$coefficients[2,]
}
res <- as.data.frame(res)
res$fdr <- p.adjust(res$p,method="fdr")
sum(res$fdr < 0.05)
sum(res$fdr < 0.1)
sum(res$fdr < 0.2)
sum(res$p < 0.05)
res <- as.data.frame(cbind(vmrs2,res))
colnames(res)[1:3] <- c("seqnames","start","end")
write.csv(res,"dmr.csv")

# GO enrich
dmr <- plyranges::as_granges(res[res$p < 0.05,1:3])
go <- great(dmr,"GO:BP","RefSeq:hg38", background=paste0("chr", 1:22))
tb  <- getEnrichmentTable(go)
write.csv(tb,"dmr_great_p0.05.csv")
