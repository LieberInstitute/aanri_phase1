"""
This script is used to plot GO analysis produced by GOATOOLS.
It communicates with R, so rpy2 needs to be installed.
"""

import numpy as np
import pandas as pd
from rpy2.robjects.conversion import localconverter, py2rpy
from rpy2.robjects import pandas2ri, r, globalenv, default_converter

def get_top_GO(tissue, direction):
    if direction == "ALL":
        label = "%s.xlsx" % tissue.lower().replace(" ", "")
    elif direction == "UP":
        label = "%s_up.xlsx" % tissue.lower().replace(" ", "")
    else:
        label = "%s_down.xlsx" % tissue.lower().replace(" ", "")
    fn = "GO_analysis_mash_%s" % (label)
    df = pd.read_excel(fn).sort_values('p_uncorrected').head(10)
    df['Log10'] = -np.log10(df['p_fdr_bh'])
    df['Tissue'] = tissue
    df["Direction"] = direction
    return df


def go_df_for_plotting():
    df = pd.DataFrame()
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        for direction in ["ALL", "UP", "DOWN"]:
            df = pd.concat([df, get_top_GO(tissue, direction)], axis=0)
    fac = []
    for ii in range(df.shape[0]):
        xx, yy = df[['ratio_in_study']].iloc[ii, 0].split('/')
        zz, tt = df[['ratio_in_pop']].iloc[ii, 0].split('/')
        fac.append((int(xx) / int(yy)) / (int(zz) / int(tt)))
    df['geneRatio'] = fac
    return df.sort_values('p_uncorrected')


def plot_go():
    pandas2ri.activate()
    globalenv['r_df'] = go_df_for_plotting()
    r("""
    library(dplyr)
    library(ggplot2)
    save_plot <- function(p, fn, w, h){
        for(ext in c('.svg', '.pdf')){
            ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    plot_GO <- function(df){
        cbPalette <- ggpubr::get_palette(palette = "npg", 4)
        gg1 = ggplot(df, aes(x=Log10, y=name, color=Tissue, size=geneRatio)) +
            geom_point(shape=18, alpha=0.8) + labs(y='', x='-Log10 (FDR)') +
            scale_colour_manual(name="Brain Region", values=cbPalette,
                                labels=c("Caudate", "Dentate Gyrus",
                                         "DLPFC", "Hippocampus")) +
            geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
            scale_size_continuous(range = c(2.5, 12)) +
            theme_bw(base_size=15) +
            theme(axis.title=element_text(face='bold'),
                  strip.text=element_text(face='bold'))
        return(gg1)
    }
    ## Plotting
    for(direction in c("ALL", "UP", "DOWN")){
        dx = r_df %>% filter(Direction == direction)
        gg = plot_GO(dx)
        config = list("DOWN"="aa_bias", "UP"="ea_bias", "ALL"="all")
        save_plot(gg, paste0("ancestry_mash_GO_top10_stacked_",
                             config[[direction]]), 8, 6)
    }
    """)


def main():
    plot_go()


if __name__ == '__main__':
    main()
