## Plotting analysis
library(dplyr)
library(ggpubr)

save_plot <- function(p, fn, w, h){
    for(ext in c('.png', '.pdf')){
        ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

#### MAIN
mds <- data.table::fread("../../../input/genotypes/mds/_m/LIBD_Brain_TopMed.mds") %>%
    select(FID, C1, C2)
pheno_df <- data.table::fread("../../../input/phenotypes/merged/_m/merged_phenotypes.csv")
dg     <- pheno_df %>% filter(Region == "DentateGyrus") %>%
    select(BrNum, Race) %>% mutate(DG="Yes") %>% distinct(.keep_all=TRUE)
others <- pheno_df %>% filter(!(BrNum %in% dg$BrNum)) %>%
    select(BrNum, Race) %>% mutate(DG="No") %>% distinct(.keep_all=TRUE)
new_pheno <- bind_rows(others, dg) %>% distinct(.keep_all=TRUE)
df  <- new_pheno %>% inner_join(mds, by=c("BrNum"="FID")) %>%
    distinct(BrNum, .keep_all=TRUE)

                                        # Plotting
pp = ggpubr::ggscatter(df, x="C1", y="C2", color="DG", size=1,
                       xlab="Component 1", ylab="Component 2",
                       ggtheme=ggpubr::theme_pubr(base_size=15))
fn = "mds_AAonly_plotting"
save_plot(pp, fn, 6, 6)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
