# compute sd of DNAm level for each CpG

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
args = commandArgs(trailingOnly=TRUE)

chr <- args[1]

# load BS object
BSobj = loadHDF5SummarizedExperiment('/dcl02/lieber/WGBS/psychENCODE_szControl/Batch3/bs_objs/batch3_combined/', prefix='CpG')
load("/dcl01/lieber/ajaffe/Kira/SczWGBS/blacklist.rda")
bb = findOverlaps(BSobj, blacklist)
BSobj = BSobj[-queryHits(bb),]
BSobj = BSobj[seqnames(BSobj) == chr,] #keep chr

# remove CT snps at CpG sites
f_snp <- paste0("/dcs04/lieber/statsgen/shizhong/AANRI/DEM2/snps_CT/",chr)
snp <- fread(f_snp,header=F,data.table=F)
idx <- is.element(start(BSobj), snp[,1])
BSobj <- BSobj[!idx,]

# keep AA only
f_ances <- "/dcs04/lieber/statsgen/shizhong/AANRI/structure/structure_CEU_AFR/structure.out_ancestry_proportion_raceDemo_compare"
f_demo <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC"
ances <- read.table(f_ances,header=T)
demo  <- read.table(f_demo,header=T)
id <- intersect(intersect(ances$id[ances$group == "AA"],colData(BSobj)$brnum),demo$BrNum[demo$Dx == "Control" & demo$Age >= 17])
BSobj <- BSobj[,is.element(colData(BSobj)$brnum,id)] 

# exlcude low coverage sites
cov=getCoverage(BSobj)
n <- length(colData(BSobj)$brnum)
keep <- which(rowSums2(cov >= 5) >= n * 0.8) 
BSobj <- BSobj[keep,]

# Compute means and sds
M=as.matrix(getMeth(BSobj))
sds <- rowSds(M)
means <- rowMeans2(M)
save(sds,means,BSobj,file=paste0("./out/",chr,".rda"))
