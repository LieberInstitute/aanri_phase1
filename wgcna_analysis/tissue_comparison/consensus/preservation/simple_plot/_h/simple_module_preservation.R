## Display module preservation results without In-Group Proportion
library(WGCNA)

load(file="../../_m/04.RData", verbose=TRUE)

## Print preservation summary statistics
ref = 2; test = 1
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1],
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1],
               mp$preservation$Z[[ref]][[test]][, -1])
## Compare preservation to quality
print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)))
## Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]
                                        # Leave grey and gold out
plotMods = !(modColors %in% c("grey", "gold"))
                                        # Text labels for points
text = modColors[plotMods]
                                        # Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2],
                 mp$preservation$Z[[ref]][[test]][,2])
                                        # Main titles for the plot
mains = c("Preservation Median Rank", "Preservation Zsummary")
##sizeGrWindow(10, 5);
pdf(file="modulePreservation_Zsummary_medianRank.pdf", width=10, height=5)
par(mfrow = c(1,2)); par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2){
    min = min(plotData[, p], na.rm = TRUE);
    max = max(plotData[, p], na.rm = TRUE);
                                        # Adjust ploting ranges appropriately
    if (p==2){
        if (min > -max/10) min = -max/10
        ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
    } else
        ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
    plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1,
         bg = modColors[plotMods], pch = 21, main = mains[p], cex = 2.4,
         ylab = mains[p], xlab = "Module size", log = "x",
         ylim = ylim, xlim = c(10, 2000), cex.lab = 1.2,
         cex.axis = 1.2, cex.main =1.4);
    labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1,
                offs = 0.08);
                                        # For Zsummary, add threshold lines
    if (p==2){
        abline(h=0); abline(h=2, col = "blue", lty = 2)
        abline(h=10, col = "darkgreen", lty = 2)
    }
}
dev.off()
