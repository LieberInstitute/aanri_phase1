#!/usr/bin/env Rscript
## @author = "Kynon JM Benjamin"

suppressMessages({
    library(dplyr)
    library(bsseq)
    library(optparse)
    library(GenomicRanges)
})

get_position <- function(chr, start, end){
    buf = 2e4 # gene position +/- 20 kb
    gr <- GRanges(seqnames=chr, ranges=IRanges(start-buf, end+buf))
    return(gr)
}

map_tissues  <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="dlpfc")[[tissue]])
}

get_bsobj <- function(tissue, chrom, methyl_dir){
    load(paste0(cpg_base,map_tissues(tissue),"_",
                chrom,"BSobj_GenotypesDosage.rda"))
    return(BSobj)
}

memBSobj <- memoise::memoise(get_bsobj)

get_overlapHits <- function(chr, start, end, tissue, methyl_dir){
    gr <- get_position(chr, start, end)
    oo = findOverlaps(memBSobj(tissue, chr, methyl_dir), gr)
    return(memBSobj(tissue, chr, methyl_dir)[queryHits(oo),])
}

save_methylation <- function(chr, start, end, outdir, tissue, methyl_dir){
    regionBS <- get_overlapHits(chr, start, end, tissue, methyl_dir)
    annot <- granges(regionBS) %>% as.data.frame %>%
        mutate(cpg_id=paste(seqnames, start, sep="_"))
    methyl <- assays(regionBS)$M %>% as.data.frame %>%
        mutate(cpg_id=annot$cpg_id) %>%
        tibble::column_to_rownames("cpg_id")
    colnames(methyl) <- colData(regionBS)$brnum
    data.table::fwrite(x=methyl,
                       file=paste0(outdir, "/cpg_methylation_data.csv"),
                       sep=',', row.names=TRUE)
}

opt <- parse_args(
    OptionParser(
        usage = "%prog [options]",
        description="Extract methylation data based on gene position.",
        option_list=list(
            make_option(c("-c","--chrom"), action="store",
                        default = "chr1", type="character",
                        help = "Chromosome name [default %default]"),
            make_option(c("-s","--start"), action="store",
                        default=NA, type="integer",
                        help = "Gene start position."),
            make_option(c("-e", "--end"), action="store",
                        default=NA, type="integer",
                        help = "Gene end position."),
            make_option(c("-m", "--methyl_dir"), action="store",
                        default=NA, type="character",
                        help = "PATH to the normalized methylation BS object."),
            make_option(c("-t", "--tissue"), action="store",
                        default="caduate", type="character",
                        help="Brain region to analyze [default %default]"),
            make_option(c("-o", "--output"), action="store",
                        default = "output", type="character",
                        help = "Output directory [default %default]")
        )
    )
)

save_methylation(opt$chrom, opt$start, opt$end, opt$output,
                 tolower(opt$tissue), opt$methyl_dir)
