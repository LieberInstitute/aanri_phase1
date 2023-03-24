library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
library(rGREAT)

setwd("out")

# read files
fs <- list.files(path = ".", pattern = "^p",recursive = T, full.names = T)
fs <- fs[!grepl("pc",fs)]
raw <-  lapply(fs, function(x) fread(x,data.table=F))
v <- rbindlist(raw)
colnames(v) <- c("chr","start","sd")

# sd cutoff
sdCut <- quantile(v[,3], prob = 0.99, na.rm = TRUE)

# get vmr
vmrs <- c()
for(i in 1:22){
	cat(i,"\n")
	v2 <- v[v$chr==paste0("chr",i),]
	v2 <- v2[order(v2$start)]
	isHigh <- rep(0, nrow(v2))
	isHigh[v2$sd > sdCut] <- 1
	vmrs0 <- bsseq:::regionFinder3(isHigh, as.character(v2$chr), v2$start, maxGap = 1000)$up
	vmrs <- rbind(vmrs,vmrs0)
}
vmr <- vmrs[vmrs$n > 5,1:3]
write.table(vmr,"vmr.bed",col.names=F,row.names=F,sep="\t",quote=F)

# GO enrich
colnames(vmr) <- c("seqnames","start","end")
vmr <- plyranges::as_granges(vmr)
res <- great(vmr,"GO:BP","RefSeq:hg38", background=paste0("chr", 1:22))
tb  <- getEnrichmentTable(res)
write.csv(tb,"vmr_great.csv")
