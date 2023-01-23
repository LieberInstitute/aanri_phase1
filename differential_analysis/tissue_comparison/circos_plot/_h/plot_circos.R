## Generate circlized plot
suppressPackageStartupMessages({
    library(dplyr)
    library(biomaRt)
    library(circlize)
    library(ComplexHeatmap)
})

get_annotation <- function(){
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    biomart = getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                 'start_position', 'end_position'),
                    mart=ensembl)
    return(biomart)
}
memMART <- memoise::memoise(get_annotation)

extract_bed <- function(tissue){
    fn = "../../summary_table/_m/BrainSeq_ancestry_4features_4regions.txt.gz"
    bed = data.table::fread(fn) %>% filter(Type == "Gene", Tissue == tissue) %>%
        mutate(ensembl_gene_id=gsub("\\..*", "", gencodeID)) %>%
        inner_join(memMART(), by=c("ensembl_gene_id")) %>%
        dplyr::select(chromosome_name,start_position,
                      end_position,posterior_mean,lfsr) %>%
        mutate(chromosome_name=paste0('chr', chromosome_name))
    bed_EA = bed %>% filter(posterior_mean > 0)
    bed_AA = bed %>% filter(posterior_mean < 0)
    return(list("EA"=bed_EA, "AA"=bed_AA))
}

plot_circos_4tissue <- function(caudate, dlpfc, hippo, gyrus){
    lgd_points = Legend(at=c("Increased EA proportion", "Increased AA proportion"),
                        type="points", legend_gp=gpar(col = c("red", "blue")),
                        title_position="topleft", title="Genetic Ancestry",
                        background="#FFFFFF")
    circos.clear() # clear plot if there is any
    circos.par("start.degree" = 0, "cell.padding" = c(0, 0, 0, 0),
               "track.height" = 0.15) # rotate 90 degrees
    # initialize with ideogram
    # use hg38, default is hg19
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicTrack(caudate, bg.border="#E64B35FF",
                        bg.col=add_transparency("#E64B35FF", transparency=0.7),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("red", "blue")[i], ...)
    })
    circos.genomicTrack(gyrus, bg.border="#4DBBD5FF",
                        bg.col=add_transparency("#4DBBD5FF", transparency=0.7),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("red", "blue")[i], ...)
    })
    circos.genomicTrack(dlpfc, bg.border="#00A087FF",
                        bg.col=add_transparency("#00A087FF", transparency=0.7),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("red", "blue")[i], ...)
    })
    circos.genomicTrack(hippo, bg.border="#3C5488FF",
                        bg.col=add_transparency("#3C5488FF", transparency=0.7),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("red", "blue")[i], ...)
    })
    draw(lgd_points, x=unit(5, "mm"), y=unit(5, "mm"), just=c("left", "bottom"))
}

####### MAIN
main <- function(){
    caudate <- extract_bed("Caudate")
    gyrus   <- extract_bed("Dentate Gyrus")
    dlpfc   <- extract_bed("DLPFC")
    hippo   <- extract_bed("Hippocampus")
                                        # plot
    pdf(file = paste0("significant_circos_plot_4regions.pdf"),
        width = 10, height = 10)
    plot_circos_4tissue(caudate, dlpfc, hippo, gyrus)
    dev.off()
}

main()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
