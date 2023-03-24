library(rGREAT) 

# load vmrs
load("../out/chr1_vmr.rda")
meth2 <- meth
vmrs2 <- vmrs
for(i in 2:22){
	load(paste0("../out/chr",i,"_vmr.rda"))
	meth2 <- rbind(meth2,meth)
	vmrs2 <- rbind(vmrs2,vmrs)
}
p <- t(meth2)
ind <- rownames(p)
colnames(p) <- paste0(vmrs2[,1],"_",vmrs2[,2],"_",vmrs2[,3])

# replace missing data by col means
for(i in 1:ncol(p)){
	p[is.na(p[,i]),i] <- mean(p[,i],na.rm=T)
}

# read local ancesty
f_ances <- "./out/local_ancesty.csv"
ances <- read.csv(f_ances,header=T)
# remove regions without ancestry estimates
idx <- colSums(is.na(ances))>=20
ances <- ances[,!idx]
ances <- ances[match(ind,ances$X),]

# read covariates
f_demo <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC"
f_sv <- "../out/pc.csv"
demo  <- read.table(f_demo,header=T)
demo <- demo[match(ind,demo$BrNum),]
sv <- read.csv(f_sv)
sv <- sv[match(ind,sv$X),-1]

# match regions
gr <- intersect(colnames(ances),colnames(p))
p <- p[,match(gr,colnames(p))]
ances <- ances[,match(gr,colnames(ances))]

# differential DNAm for Afr proportion
dmr <- matrix(0,nrow=ncol(p),ncol=4)
residual <- p
residual[] <- NA
for(i in 1:ncol(p)){
	if(! i %% 100){
		cat(i,"\n")
	}
	# regress out local ancestry, age, sex, top 5 PC
	d <- as.data.frame(cbind(y=p[,i],ances[,i]/2,sv[,1:5],demo$Age,as.factor(demo$Sex)))
	model = lm(y ~ .,data=d)
	dmr[i,] = summary(model)$coefficients[2,]
	residual[,i] <- resid(model)
}
colnames(dmr) <- c("beta","se", "t", "p")
dmr <- as.data.frame(dmr)
dmr$fdr <- p.adjust(dmr$p,method="fdr")
sum(dmr$fdr<0.05)
sum(dmr$fdr<0.1)
sum(dmr$p<0.05)

# output
chr <- sapply(colnames(p),function(x){strsplit(x,"_")[[1]][1]})
start <- sapply(colnames(p),function(x){strsplit(x,"_")[[1]][2]})
end <- sapply(colnames(p),function(x){strsplit(x,"_")[[1]][3]})
bed <- data.frame(seqnames=chr,start=as.numeric(start),end=as.numeric(end))
res <- cbind(bed,dmr)
write.table(res,"./out/dmr.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(residual,"./out/vmr_residual.txt",col.names=T,row.names=T,sep="\t",quote=F)

# monitor local and global ancestry for each DMR region
f_ances <- "/dcs04/lieber/statsgen/shizhong/AANRI/structure/structure_CEU_AFR/structure.out_ancestry_proportion_raceDemo_compare"
ancesg <- read.table(f_ances,header=T)
ancesg <- ancesg[match(ind,ancesg$id),]
dmr_global <- read.csv("../out/dmr.csv")
dmr_global$region <- paste(dmr_global[,2],dmr_global[,3],dmr_global[,4],sep="_")
pdf("./out/DMR_global_local_compare_by_region.pdf")
par(mfrow=c(2,2))
dmr_local <- res[res$fdr<0.05,]
for(i in 1:nrow(dmr_local)){
	idx <- which(gr == rownames(dmr_local)[i])
	idx2 <- which(dmr_global$region == rownames(dmr_local)[i])
	main <- paste(rownames(dmr_local)[i],
	    # beta, se and t
		paste0("local: ", round(res[idx,4],2),",",round(res[idx,5],2),",",round(res[idx,6],2)),
		paste0("global: ", round(dmr_global[idx2,6],2),",",round(dmr_global[idx2,7],2),",",round(dmr_global[idx2,8],2)),sep="\n")
	plot(ances[,idx]/2,p[,idx],xlim=c(0,1),col="red",xlab="ancestry",ylab="DNAm",main=main)
	points(ancesg$Afr,p[,idx],col="blue")
}
dev.off()	

# global effect size compare
pdf("./out/DMR_global_local_compare_overall.pdf")
region <- intersect(dmr_global$region, gr)
dmr_global2 <- dmr_global[match(region,dmr_global$region),]
res2 <- res[match(region,gr),]
cor(res2$beta,dmr_global2$beta)
plot(res2$beta,dmr_global2$beta,xlab="Local ancestry effects",ylab="Global ancestry effects")
dev.off()

# GO analysis
dmr <- plyranges::as_granges(res[res$fdr < 0.05,1:3])
go <- great(dmr,"GO:BP","RefSeq:hg38", background=paste0("chr", 1:22))
tb  <- getEnrichmentTable(go)
write.csv(tb,"./out/dmr_great_fdr0.05.csv")