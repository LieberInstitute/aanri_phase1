##### Review correlation between VMR PCs and SNP PCs

library(dplyr)

get_vmr_pcs <- function(region){
    region_lt <- list("Caudate"="caud",
                      "DLPFC"="dlpfc",
                      "Hippocampus"="hippo")
    fn <- paste0("/dcs04/lieber/statsgen/shizhong/",
                 "AANRI/VMR3/99/", region_lt[[region]],
                 "/aa2/out/pc.csv")
    return(data.table::fread(fn) |>
           select(V1, PC1, PC2, PC3, PC4, PC5) |>
           tidyr::pivot_longer(!V1, names_to="VMR_PCs",
                               values_to="VMR_var") |>
           as.data.frame() |> mutate(Region=region))
}

merge_vmrs <- function(){
    dt = list()
    for(region in c("Caudate", "DLPFC", "Hippocampus")){
        dt[[region]] <- get_vmr_pcs(region)
    }
    return(bind_rows(dt))
}

get_snp_pcs <- function(){
    fn <- here::here("input/genotypes/mds/_m/TOPMed_LIBD_AA.mds")
    return(data.table::fread(fn) |> select(-IID, -SOL) |>
           rename_at(.vars = vars(starts_with("C")),
                     function(x){sub("C", "snpPC", x)}) |>
           tidyr::pivot_longer(!FID, names_to="SNP_PCs",
                               values_to="SNP_var") |>
           as.data.frame())
}

merge_data <- function(){
    return(inner_join(get_snp_pcs(), merge_vmrs(),
                      by=c("FID"="V1"),
                      relationship="many-to-many") |>
           rename("BrNum"="FID"))
}

#### Main
for(region in c("Caudate", "DLPFC", "Hippocampus")){
    print(region)
    est_fit <- merge_data() |> filter(Region == region) |>
        group_by(SNP_PCs, VMR_PCs) |>
        do(fitEST = broom::tidy(lm(SNP_var ~ VMR_var, data = .))) |>
        tidyr::unnest(fitEST) |> filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf, breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf+10**(-100)))

    print(est_fit |> count(p.bonf.cat))
    print(est_fit |> filter(p.value < 0.05))
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

