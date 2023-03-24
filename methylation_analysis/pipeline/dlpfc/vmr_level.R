library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library(scales)
args = commandArgs(trailingOnly=TRUE)

chr <- args[1]

setwd("out")

# read vmr
vmr <- read.table("vmr.bed",header=F)
vmrs <- vmr[vmr[,1] ==  chr,]

# load BSobj
load(paste0(chr,".rda"))

# compute DNAm levels for each vmr
meth <- matrix(NA,nrow=nrow(vmrs),ncol=nrow(colData(BSobj)))
n <- rep(0,nrow(vmrs)) # nb of CpG within each vmr
for(i in 1:nrow(vmrs)){
	cat(i,"\n")
	idx <- start(BSobj) >= vmrs[i,2] & start(BSobj) <= vmrs[i,3]
	n[i] <- sum(idx)
	BSregion <- BSobj[idx,]
	m <- getCoverage(BSregion,type="M")
	cov <- getCoverage(BSregion,type="Cov")
	meth[i,] <- colSums(m)/colSums(cov)
}
colnames(meth) <- colData(BSobj)$brnum
vmrs <- cbind(vmrs,n)
save(vmrs,meth,file=paste0(chr,"_vmr.rda"))

