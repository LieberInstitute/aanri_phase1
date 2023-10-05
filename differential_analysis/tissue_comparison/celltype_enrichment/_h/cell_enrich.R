fdr_cutoff <- 0.05

# read deg results
infile <- here::here("differential_analysis/tissue_comparison",
                     "summary_table/_m",
                     "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
dat <- read.table(gzfile(infile),header=T,sep="\t")

# keep only genes
dat <- dat[dat$Type=="Gene",]
dat$Tissue[dat$Tissue=="Dentate Gyrus"] <- "Dentate_Gyrus"

# background genes
bg <- unique(dat$Symbol[dat$Type=="Gene"])
bg <- bg[!is.na(bg)]

# get degs
tissue <- unique(dat$Tissue)
geneList <- list()
k=0
for(i in 1:length(tissue)){
	k <- k +1
	idx <- dat$Tissue == tissue[i] & dat$lfsr <= fdr_cutoff
	geneList[[k]] <- unique(dat$Symbol[idx])
	idx <- dat$Tissue == tissue[i] & dat$posterior_mean > 0  & dat$lfsr <= fdr_cutoff
	geneList[[k + 1]] <- unique(dat$Symbol[idx])
	idx <- dat$Tissue == tissue[i] & dat$posterior_mean < 0  & dat$lfsr <= fdr_cutoff
	geneList[[k + 2]] <- unique(dat$Symbol[idx])
	names(geneList)[k:(k+2)] <- paste0(tissue[i],c("","_up","_down"))
	k <- k + 2
}

# overlap cell type genes
dat <- read.table("/dcs04/lieber/statsgen/shizhong/AANRI/ldsc3/celltype/500kb_strand/bedfiles/Zeisel_single_cell",header=T,sep="\t")
cell <- colnames(dat)[-1]
res <- c()
for(i in 1:length(cell)){
	cat(i,"\n")
	idx <- which(colnames(dat) == cell[i])
	genes.cell <- dat$geneName[dat[,idx]==1]
	for(k in 1:length(geneList)){
		genes.deg <- geneList[[k]]
		genes.bg <- bg[!is.element(bg,genes.deg)]
		n11 <- length(intersect(genes.deg,genes.cell))
		n12 <- length(genes.deg) - n11
		n01 <- length(intersect(genes.bg,genes.cell))
		n02 <- length(genes.bg) - n01
		mat <- matrix(c(n11,n12,n01,n02),2,2)
		ft <- fisher.test(mat)		
		res <- rbind(res,c(cell[i],names(geneList[k]),c(n11,n12,n01,n02),ft$estimate,ft$p.value))
	}	
}
res <- as.data.frame(res)
write.csv(res,"cell_enrich.csv")
