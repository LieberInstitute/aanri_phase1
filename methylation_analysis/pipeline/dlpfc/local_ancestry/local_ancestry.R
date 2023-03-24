library(data.table)
args = commandArgs(trailingOnly=TRUE)

# local ancestry of the gene
dat <- fread(args[1],data.table=F)
dat <- dat[,-(1:4)]

# index for AFR haplotypes
idx <- seq_len(ncol(dat))%%2
# mean for each AFR haplotype
hap <- colMeans(dat[,idx==1])
# sum over two AFR haplotypes
idx <- seq_len(length(hap))%%2
hap2 <- hap[idx==1] + hap[idx==0]

# sample id
id <- gsub(":::hap1:::AFR","",names(hap[idx==1]))

res <- data.frame(ID=id,AFR=hap2)	
write.table(res,args[2],col.names=T,row.names=F,sep="\t",quote=F)


