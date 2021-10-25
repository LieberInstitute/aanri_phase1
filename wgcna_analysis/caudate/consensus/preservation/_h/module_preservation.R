## Examine module preservation
suppressPackageStartupMessages({
    library(WGCNA)
    library(clusterRepro)
})

## Load networks
load("../../_m/01.RData", verbose=TRUE)
load("../../_m/02.RData", verbose=TRUE)
load("../../_m/03.RData", verbose=TRUE)

## Generate color list
colorAA = moduleColors
colorEA = moduleColors
colorList = list(colorAA, colorEA)
names(colorList) = shortLabels

## Module preservation
enableWGCNAThreads()
mp = modulePreservation(multiExpr, colorList,
                        parallelCalculation=FALSE,
                        referenceNetworks = c(1:2),
                        loadPermutedStatistics = FALSE,
                        corFnc="bicor", networkType="signed",
                        nPermutations=200, verbose=3, randomSeed=13)
save(mp, file="04.RData")

## Print preservation summary statistics
ref = 2
test = 1
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1],
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1],
               mp$preservation$Z[[ref]][[test]][, -1])
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

## Calculate In-Group Proportion
impExpr = list()
nSets = exprSize$nSets
for(set in 1:nSets){
    impExpr[[set]] = list(data=t(impute::impute.knn(t(multiExpr[[set]]$data))$data))
}
eigengenes = list()
for(set in 1:nSets){
    eigengenes[[set]] = multiSetMEs(impExpr, universalColors=colorList[[set]],
                                    excludeGrey=TRUE)
    for(ss in 1:nSets){
        rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExpr[[ss]]$data)
    }
}
                                        # In-Group Proportion
cr = list()
set.seed(13)
for(ref in 1:nSets){
    cr[[ref]] = list()
    for(test in 1:nSets){
        cr[[ref]][[test]] = clusterRepro(Centroids=as.matrix(eigengenes[[ref]][[test]]$data),
                                         New.data=as.matrix(impExpr[[test]]$data),
                                         Number.of.permutations=1000)
        collectGarbage()
    }
}
save(cr, file="05.RData")

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
