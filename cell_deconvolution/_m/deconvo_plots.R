
library(SummarizedExperiment)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(broom)
library(viridis)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(patchwork)
library(jaffelab)


## Load Bisque results & res data
load("est_prop_Bisque.v2.Rdata",verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution/data/cell_colors.Rdata", verbose = TRUE)

pd_fn <- list(dg = "dg_phenotypes.csv", 
              caudate = "caudate_phenotypes.csv",  
              hippo = "hippo_phenotypes.csv", 
              dlpfc = "dlpfc_phenotypes.csv")

pd_path <- "/dcl01/beer/lieber/kynon_temp/control_ancestry/input/phenotypes/_m/"
pd <- map(pd_fn, ~read.csv(paste0(pd_path, .x), row.names = "X"))
pd <- do.call("rbind", pd)

pd2 <- pd[,c("RNum", "BrNum", "Region", "Sex", "Race")]
rownames(pd2) <- NULL

long_prop <- do.call("rbind", map(est_prop_bisque, "Est.prop.long")) %>%
  separate(sample, into = c("RNum")) %>%
  left_join(pd2) %>%
  mutate(`Region + Race` = paste(Region, Race))
  
long_prop %>% filter(prop > 0.99)
#     RNum cell_type prop  BrNum  Region Sex Race Region + Race
# 1 R13321     Oligo    1 Br1492 Caudate   M   AA    Caudate AA
# 2  R4965     Oligo    1 Br1804   HIPPO   M   AA      HIPPO AA
# 3  R5460     Oligo    1 Br1455   HIPPO   M   AA      HIPPO AA
# 4  R3464     Oligo    1 Br1811   DLPFC   F   AA      DLPFC AA
# 5  R4222     Inhib    1 Br1647   DLPFC   M CAUC    DLPFC CAUC

## Boxplots
boxplot_region <- long_prop %>%
  ggplot(aes(x = cell_type, y = prop, fill = Region)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill = 'Brain Region') +
  theme_bw(base_size = 15)

ggsave(boxplot_region, filename = "plots/cellType_boxplot_region.png", width = 12)
ggsave(boxplot_region, filename = "plots/cellType_boxplot_region.pdf", width = 12)

boxplot_race <- long_prop %>%
  ggplot(aes(x = Region, y = prop, fill = Race)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Primary Dx') +
  theme_bw(base_size = 15) +
  facet_wrap(~cell_type, scale = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(boxplot_race, filename = "plots/cellType_boxplot_race.png", width = 12, height = 12)
ggsave(boxplot_race, filename = "plots/cellType_boxplot_race.pdf", width = 10)

#### composition barplot ####
comp_barplot <- plot_composition_bar(long_prop, x_col = "Region") +
  scale_fill_manual(values = cell_colors)+
  theme_bw(base_size = 15)

ggsave(comp_barplot, filename = here("plots","bisque_composition_barplot.png"))
ggsave(comp_barplot, filename = here("plots","bisque_composition_barplot.pdf"))

long_prop2 <- long_prop %>%
  mutate(`Region + Dx` = paste(BrainRegion, PrimaryDx))

long_prop2 %>% count(`Region + Dx`)

comp_barplot_race <- plot_composition_bar(long_prop, x_col = "Region + Race") +
  scale_fill_manual(values = cell_colors)+
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(comp_barplot_race, filename = here("plots","bisque_composition_barplot_race.png"), height = 8, width = 12)
ggsave(comp_barplot_race, filename = here("plots","bisque_composition_barplot_race.pdf"), height = 8, width = 12)

#### Cor with qSV ####
qSV_path = "/dcl01/beer/lieber/kynon_temp/control_ancestry/input/phenotypes/_m/"
qSV_fn <- list(dg = "qSV_dg.csv", 
              caudate = "qSV_caudate.csv",  
              hippo = "qSV_hippo.csv", 
              dlpfc = "qSV_dlpfc.csv")

qSV_mat <- map(qSV_fn, ~read.csv(paste0(qSV_path, .x), row.names = "X"))
map(qSV_mat, dim)

qSV_long <- map(qSV_mat, ~as.matrix(.x) %>% melt()  %>% rename(sample = Var1, qSV = Var2, qSV_value = value))

## Bind with qSV table

est_prop_bisque$dg$Est.prop.long$sample <- ss(as.character(est_prop_bisque$dg$Est.prop.long$sample), "_")
  
est_prop_qsv <- map2(est_prop_bisque, qSV_long, ~left_join(.x$Est.prop.long, .y, by = "sample"))
map(est_prop_qsv, head)

#### Calculate p-values ####
prop_qSV_fit <- map(est_prop_qsv, ~.x %>% 
                      group_by(cell_type, qSV) %>%
                      do(fitQSV = tidy(lm(prop ~ qSV_value, data = .))) %>%
                      unnest(fitQSV) %>%
                      filter(term == "qSV_value") %>%
                      mutate(p.bonf = p.adjust(p.value, "bonf"),
                             p.bonf.sig = p.bonf < 0.05,
                             p.bonf.cat = cut(p.bonf,
                                              breaks = c(1,0.05, 0.01, 0.005, 0),
                                              labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05")
                             ),
                             p.fdr = p.adjust(p.value, "fdr"),
                             log.p.bonf = -log10(p.bonf)))

map(prop_qSV_fit, ~.x %>% count(p.bonf.cat))
# $dg
# # A tibble: 3 × 2
# p.bonf.cat     n
# <fct>      <int>
#   1 <= 0.005      24
# 2 <= 0.05        3
# 3 > 0.05        33
# 
# $caudate
# # A tibble: 4 × 2
# p.bonf.cat     n
# <fct>      <int>
#   1 <= 0.005      39
# 2 <= 0.01        1
# 3 <= 0.05        6
# 4 > 0.05        84
# 
# $hippo
# # A tibble: 4 × 2
# p.bonf.cat     n
# <fct>      <int>
#   1 <= 0.005      53
# 2 <= 0.01        1
# 3 <= 0.05        6
# 4 > 0.05        80
# 
# $dlpfc
# # A tibble: 4 × 2
# p.bonf.cat     n
# <fct>      <int>
#   1 <= 0.005      42
# 2 <= 0.01        1
# 3 <= 0.05        5
# 4 > 0.05        42

walk2(prop_qSV_fit, names(prop_qSV_fit), ~write_csv(.x, paste0("qSV_vs_prop_fit-",.y,".csv")))

#### Tile plots ####
walk2(prop_qSV_fit, names(prop_qSV_fit), function(fit, name){
  
  tile_plot_val <- fit %>%
    ggplot(aes(x = qSV, y = cell_type, fill = log.p.bonf)) +
    geom_tile(color = "grey") +
    geom_text(aes(label = ifelse(p.bonf.sig, format(round(-log10(p.bonf),1), nsmall = 1),""),
                  color = log10(p.bonf)), size = 3, fontface = "bold",
              show.legend = F)+
    # scale_color_viridis(option = "magma") +
    scale_fill_viridis(name = "-log10(p-value Bonf)", option = "magma", direction = -1) +
    scale_y_discrete(limits = rev) +
    labs(title = paste(name, "p-values cell-type prop~qSV"), x = 'Cell Type', color = "p-value Bonf\nsignificance") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(plot = tile_plot_val,
         filename = here("plots",paste0("qSV_prop_fit_tileVal-",name,".png")),
         width = 10)
  
  ggsave(plot = tile_plot_val,
         filename = here("plots",paste0("qSV_prop_fit_tileVal-",name,".pdf")),
         width = 10)
})


#### Create scatter plots ####
# sig_colors2 <- c(brewer.pal(3, "Set1"),"black")
# sig_colors2 <- c("#440154","#31688E", "#35B779","black")
# names(sig_colors2) <- levels(prop_qSV_fit$p.bonf.cat)
# 
# est_prop_qsv_fit <- left_join(est_prop_qsv, prop_qSV_fit)
# 
# scatter_plot <- map(regions, ~  est_prop_qsv_fit %>%
#                       filter(BrainRegion == .x) %>%
#                       ggplot(aes(x = qSV_value, y = prop, color = p.bonf.cat))+
#                       geom_point(size = .4, alpha = .2) +
#                       facet_grid(cell_type~qSV, scales = "free")+
#                       theme_bw(base_size = 10)+
#                       scale_color_manual(values = sig_colors2) +
#                       theme(legend.text = element_text(size = 15)) +
#                       guides(color = guide_legend(override.aes = list(size=5))) +
#                       labs(title = .x)
# )
# 
# 
# walk2(scatter_plot, names(scatter_plot), ~ggsave(filename = here( "plots", paste0("qSV_cellType_scatter-",.y,".png")),
#        plot = .x, width = 26, height = 10))



# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.0 Patched (2021-05-18 r80330)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2021-06-21
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.2.1    2020-12-09 [2] CRAN (R 4.1.0)
# beachmat               2.8.0    2021-05-19 [2] Bioconductor
# Biobase              * 2.52.0   2021-05-19 [2] Bioconductor
# BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor
# BiocNeighbors          1.10.0   2021-05-19 [1] Bioconductor
# BiocParallel           1.26.0   2021-05-19 [2] Bioconductor
# BiocSingular           1.8.1    2021-06-08 [1] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# bluster                1.2.1    2021-05-27 [1] Bioconductor
# broom                * 0.7.7    2021-06-13 [2] CRAN (R 4.1.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    2.5.0    2021-04-26 [2] CRAN (R 4.1.0)
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.0)
# codetools              0.2-18   2020-11-04 [2] CRAN (R 4.1.0)
# colorout             * 1.2-2    2021-05-27 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.0-1    2021-05-04 [2] CRAN (R 4.1.0)
# crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DeconvoBuddies       * 0.99.0   2021-06-14 [1] Github (lahuuki/DeconvoBuddies@d889340)
# DelayedArray           0.18.0   2021-05-19 [2] Bioconductor
# DelayedMatrixStats     1.14.0   2021-05-19 [2] Bioconductor
# digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
# dplyr                * 1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.0)
# edgeR                  3.34.0   2021-05-19 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
# GenomeInfoDb         * 1.28.0   2021-05-19 [2] Bioconductor
# GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor
# GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor
# ggplot2              * 3.3.4    2021-06-16 [2] CRAN (R 4.1.0)
# glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.1    2021-04-23 [2] CRAN (R 4.1.0)
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
# hms                    1.1.0    2021-05-17 [2] CRAN (R 4.1.0)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# igraph                 1.2.6    2020-10-06 [2] CRAN (R 4.1.0)
# IRanges              * 2.26.0   2021-05-19 [2] Bioconductor
# irlba                  2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
# jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)
# lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
# limma                  3.48.0   2021-05-19 [2] Bioconductor
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# lubridate              1.7.10   2021-02-26 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
# MatrixGenerics       * 1.4.0    2021-05-19 [2] Bioconductor
# matrixStats          * 0.59.0   2021-06-01 [2] CRAN (R 4.1.0)
# metapod                1.0.0    2021-05-19 [1] Bioconductor
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# patchwork            * 1.1.1    2020-12-17 [1] CRAN (R 4.1.0)
# pillar                 1.6.1    2021-05-16 [2] CRAN (R 4.1.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# pryr                   0.1.4    2018-02-18 [2] CRAN (R 4.1.0)
# ps                     1.6.0    2021-02-28 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
# rafalib                1.0.0    2015-08-09 [1] CRAN (R 4.1.0)
# ragg                   1.1.3    2021-06-09 [1] CRAN (R 4.1.0)
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.6    2021-01-15 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.1.0)
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
# reprex                 2.0.0    2021-04-02 [2] CRAN (R 4.1.0)
# rlang                * 0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rsvd                   1.0.5    2021-04-16 [1] CRAN (R 4.1.0)
# rvest                  1.0.0    2021-03-09 [2] CRAN (R 4.1.0)
# S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor
# ScaledMatrix           1.0.0    2021-05-19 [1] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scran                  1.20.1   2021-05-24 [1] Bioconductor
# scuttle                1.2.0    2021-05-19 [1] Bioconductor
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
# sgejobs                0.99.1   2021-05-27 [1] Github (LieberInstitute/sgejobs@f5ab0ca)
# SingleCellExperiment   1.14.1   2021-05-21 [2] Bioconductor
# sparseMatrixStats      1.4.0    2021-05-19 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# stringi                1.6.2    2021-05-17 [2] CRAN (R 4.1.0)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor
# systemfonts            1.0.2    2021-05-11 [2] CRAN (R 4.1.0)
# textshaping            0.3.5    2021-06-09 [1] CRAN (R 4.1.0)
# tibble               * 3.1.2    2021-05-16 [2] CRAN (R 4.1.0)
# tidyr                * 1.1.3    2021-03-03 [2] CRAN (R 4.1.0)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# utf8                   1.2.1    2021-03-12 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# viridis              * 0.6.1    2021-05-11 [2] CRAN (R 4.1.0)
# viridisLite          * 0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
# XVector                0.32.0   2021-05-19 [2] Bioconductor
# zlibbioc               1.38.0   2021-05-19 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
