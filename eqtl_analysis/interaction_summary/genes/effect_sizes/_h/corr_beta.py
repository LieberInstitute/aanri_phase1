"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import session_info
import pandas as pd
from functools import lru_cache
from scipy.stats import spearmanr
from rpy2.robjects.packages import importr
from rpy2.robjects import r, globalenv, pandas2ri

@lru_cache()
def get_mash_es():
    # get effect size
    fn = "../../_m/posterior_mean_allpairs_ancestry.txt.gz"
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_sig():
    fn = "../../../feature_summary/_m/" +\
        "BrainSeq_ancestry_4features_4regions.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    df["effect"] = df.gene_id + "_" +  df.variant_id
    return df[(df["Feature"] == "Gene")].copy()


@lru_cache()
def prepare_data(tissue1, tissue2):
    df1 = get_sig()[(get_sig()["Tissue"] == tissue1.replace("_", " "))].copy()
    df2 = get_sig()[(get_sig()["Tissue"] == tissue2.replace("_", " "))].copy()
    sig_genes = set(df1.effect) | set(df2.effect)
    return get_mash_es().loc[sig_genes, [tissue1, tissue2]]


@lru_cache()
def merge_data(tissue1, tissue2):
    df1 = get_sig()[(get_sig()["Tissue"] == tissue1.replace("_", " "))].copy()
    df2 = get_sig()[(get_sig()["Tissue"] == tissue2.replace("_", " "))].copy()
    sig_genes = set(df1.effect) & set(df2.effect)
    return get_mash_es().loc[sig_genes, [tissue1, tissue2]]


def corr_beta(tissue1, tissue2):
    df = prepare_data(tissue1, tissue2)
    return spearmanr(df[tissue1], df[tissue2])


def plotNsave_corr(tissue1, tissue2):
    pandas2ri.activate()
    globalenv['df'] = prepare_data(tissue1, tissue2)
    globalenv['tissue1'] = tissue1.replace("_", " ")
    globalenv['tissue2'] = tissue2.replace("_", " ")
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.png', '.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = paste0("Effect Size\n(", tissue1, ")")
    ylab = paste0("Effect Size\n(", tissue2, ")")
    tissue1 = gsub(" ", "_", tissue1); tissue2 = gsub(" ", "_", tissue2);
    fn = paste("effectsize_scatter", tissue1, tissue2, sep="_")
    pp = ggpubr::ggscatter(df, x=tissue1, y=tissue2, add="reg.line", size=1,
                           xlab=xlab, ylab=ylab,panel.labs.font=list(face="bold"),
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15), ncol=4)
    save_plot(pp, fn, 6, 6)
    ''')


def main():
    ## Calculate rho
    with open("rho_statistics.log", "w") as f:
        for tissue1 in ["Caudate", "Dentate_Gyrus", "DLPFC", "Hippocampus"]:
            for tissue2 in ["Dentate_Gyrus", "DLPFC", "Hippocampus"]:
                if tissue1 != tissue2:
                    ## Generate figure
                    plotNsave_corr(tissue1, tissue2)
                    ## Correlated effect sizes
                    rho, pval = corr_beta(tissue1, tissue2)
                    print(f"{tissue1} vs {tissue2}:\t rho = {rho:.3f}, "+\
                          f"r2 = {rho**2:.3f}, p-value = {pval:.1e}", file=f)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
