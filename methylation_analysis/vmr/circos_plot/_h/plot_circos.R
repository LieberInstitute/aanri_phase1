## Generate circlized plot
suppressPackageStartupMessages({
    library(dplyr)
    library(circlize)
    library(ComplexHeatmap)
})

extract_bed <- function(tissue, LOCAL){
    if(LOCAL){
        fn <- paste0("../../_m/aaOnly/", tissue, "_dmr_local.csv")
    } else {
        fn <- paste0("../../_m/aaOnly/", tissue, "_dmr.csv")
    }
    bed <- data.table::fread(fn) |> filter(fdr < 0.05) |>
        dplyr::rename("chr"="seqnames") |>
        dplyr::select(chr, start, end, beta)
    bed_hyper <- bed |> filter(beta > 0) |> dplyr::select(-beta)
    bed_hypo  <- bed |> filter(beta < 0) |> dplyr::select(-beta)
    return(list("hyper"=bed_hyper, "hypo"=bed_hypo))
}

plot_circos_dmr <- function(LOCAL){
                                        # Extract DMR into BED format
    caudate <- extract_bed("caudate", LOCAL)
    dlpfc   <- extract_bed("dlpfc", LOCAL)
    hippo   <- extract_bed("hippocampus", LOCAL)
                                        # Generate legend key
    lgd_points <- Legend(at=c("Hypermethylation","Hypomethylation"),
                        type="points", legend_gp=gpar(col = c("red", "blue")),
                        title_position="topleft", title="Methylation Status",
                        background="#FFFFFF")
    # clear plot if there is any
    circos.clear()
    circos.par("cell.padding" = c(0, 0, 0, 0), "track.height" = 0.15)
                                        # initialize with ideogram
                                        # use hg38, default is hg19
    circos.initializeWithIdeogram(chromosome.index=paste0("chr", 1:22),
                                  species="hg38")
    
    circos.genomicRainfall(caudate, bg.border="#E64B35FF",
                           pch=16, cex=0.6, col=c("#FF000080", "#0000FF80"),
                           bg.col=add_transparency("#E64B35FF", transparency=0.7))

    circos.genomicRainfall(dlpfc, bg.border="#4DBBD5FF",
                           bg.col=add_transparency("#4DBBD5FF", transparency=0.7),
                           pch=16, cex=0.6, col=c("#FF000080", "#0000FF80"))

    circos.genomicRainfall(hippo, bg.border="#00A087FF",
                           bg.col=add_transparency("#00A087FF", transparency=0.7),
                           pch=16, cex=0.6, col=c("#FF000080", "#0000FF80"))
    
    draw(lgd_points, x=unit(5, "mm"), y=unit(5, "mm"),
         just=c("left", "bottom"))
}

####### MAIN
                                        # Plot DMRs
pdf(file = "circos_plot_dmr_global.pdf")
plot_circos_dmr(FALSE)
dev.off()

pdf(file = "circos_plot_dmr_local.pdf")
plot_circos_dmr(TRUE)
dev.off()


#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
