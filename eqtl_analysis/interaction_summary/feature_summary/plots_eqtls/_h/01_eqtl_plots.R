## Generating eQTL boxplots
suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(ggpubr)
    library(argparse)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

plot_eqtl <- function(variant_id, gene_id, gene_name, eqtl_annot){
    y0  <- min(big_dat[, gene_id]) - 0.2
    y1  <- max(big_dat[, gene_id]) + 0.2
    titlename <- paste(gene_name, gene_id, eqtl_annot, sep="\n")
    bxp <- big_dat %>% select(all_of(c(gene_id, variant_id, "Race"))) %>%
        tidyr::drop_na() %>%
        ggboxplot(x=variant_id, y=gene_id, color="Race", add="jitter",
                  ylab="Residualized Expression", outlier.shape=NA,
                  add.params=list(alpha=0.5), alpha=0.4, legend="bottom",
                  palette="npg", ylim=c(y0,y1),
                  ggtheme=theme_pubr(base_size=20, border=TRUE)) +
        font("xy.title", face="bold") + ggtitle(titlename) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))
    return(bxp)
}

continous_plot_eqtl <- function(variant_id, gene_id, gene_name, eqtl_annot){
    y0  <- min(big_dat[, gene_id]) - 0.2
    y1  <- max(big_dat[, gene_id]) + 0.2
    titlename <- paste(gene_name, gene_id, eqtl_annot, sep="\n")
    afile <- here("input/ancestry_structure",
                  "structure.out_ancestry_proportion_raceDemo_compare")
    ancestry <- data.table::fread(afile)
    sca <- big_dat %>% select(all_of(c("BrNum", gene_id, variant_id, "Race"))) %>%
        tidyr::drop_na() %>% filter(Race == "AA") %>%
        inner_join(ancestry, by=c("BrNum"="id")) %>%
        mutate_at(variant_id, as.factor) %>%
        ggscatter(x="Afr", y=gene_id, color=variant_id, add="reg.line",
                  xlab="AA Proportion", ylab="Residualized Expression",
                  alpha=0.4, legend="bottom", ylim=c(y0,y1),
                  ggtheme=theme_pubr(base_size=20, border=TRUE)) +
        font("xy.title", face="bold") + ggtitle(titlename) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))
    return(sca)
}

#### MAIN analysis
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-t", "--tissue", type="character", default="caudate",
                    help="Tissue to be analyzed [default: %default]")
args   <- parser$parse_args()
                                        # Run analysis
feature <- "genes"; tissue  <- args$tissue

                                        # Read residualized expression
outdir  <- paste(tissue, "out", sep="/")
pfile   <- paste0(outdir, "/p")
p       <- read.table(pfile, comment.char="", header=TRUE)

                                        # gene meta
efile <- paste0(outdir, "/eqtls_top25.tsv")
e     <- read.table(efile, comment.char="", sep="\t", header=TRUE)

variants <- e[,"variant_id"]

vfile <- paste0(pfile, ".variants")
write.table(data.frame(variants), file=vfile, quote=FALSE, 
            col.names=FALSE, row.names=FALSE, sep=",")

                                        # get genotypes
plinkfile <- here("input/genotypes/all_samples/_m/TOPMed_LIBD_AA_EA")
outfile   <- paste0(pfile, "_temp")
command   <- paste("plink2 --pfile", plinkfile, "--extract", vfile,
                   "--export A", "--out", outfile)
system(command)

                                        # read sample phenotypes
sfile <- here('input/phenotypes/merged/_m/merged_phenotypes.csv')
s     <- read.table(sfile, comment.char="", header=TRUE, sep=",",
                    row.names=1)

                                        # keep samples with genotypes
p    <- p %>% tibble::column_to_rownames("gene_id") %>% t
pfam <- here("input/genotypes/all_samples/_m/TOPMed_LIBD_AA_EA.psam")
ind  <- read.table(pfam, comment.char="", header=TRUE)
id   <- intersect(rownames(p), ind[,2])
p    <- p[is.element(rownames(p), id),]

                                        # align samples phenotypes
idx  <- intersect(s$BrNum, rownames(p))
p    <- p[match(idx, rownames(p)),]
s    <- s[match(idx, s$BrNum),]
geno <- read.table(paste0(outfile,".raw"), header=TRUE)
geno <- geno[match(idx, geno$IID),] %>%
    select(-c("FID", "PAT", "MAT", "SEX", "PHENOTYPE"))
names(geno) <- gsub("_[[:upper:]]+$", "", names(geno))

                                        # prep data
big_dat <- merge(s[,c("BrNum", "Age", "Sex", "Race")], p,
                 by.x="BrNum", by.y=0) %>%
    inner_join(geno, by=c("BrNum"="IID"))

                                        # plot eQTL
plotdir <- paste(tissue, "plots", sep="/")
dir.create(plotdir)

for(j in 1:dim(e)[1]){
    cat(j,"\n")        
    eqtl       <- e[j,]; gene_id <- eqtl$gene_id
    variant_id <- eqtl$variant_id; gene_name <- eqtl$Symbol
    eqtl_annot <- paste("eQTL lfsr:", signif(eqtl$lfsr, 2))
                                        # Boxplot
    plotfile   <- paste0(plotdir, "/eqtl_", j)
    bxp        <- plot_eqtl(variant_id, gene_id, gene_name, eqtl_annot)
    save_ggplots(plotfile, bxp, 7, 7)
                                        # Scatter
    ## plotfile2  <- paste0(plotdir, "/eqtl_aa_scatter_", j)
    ## sca        <- continous_plot_eqtl(variant_id, gene_id, gene_name,
    ##                                   eqtl_annot)
    ## save_ggplots(plotfile2, sca, 7, 7)
}

                                        # clean up
outfile <- paste0(pfile, "_temp.raw")
system(paste("rm", outfile, sep=" "))
outfile <- paste0(pfile, "_temp.log")
system(paste("rm", outfile, sep=" "))

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
