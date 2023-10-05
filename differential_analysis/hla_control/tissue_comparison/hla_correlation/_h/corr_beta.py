"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import argparse
import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import spearmanr
from rpy2.robjects.packages import importr
from rpy2.robjects import r, globalenv, pandas2ri

def map_feature(feature):
    return {"genes": "Gene", "transcripts": "Transcript",
            "exons": "Exon", "junctions": "Junction"}[feature]


@lru_cache()
def get_mash_hla(feature):
    # get effect size
    fn = f"../../_m/{feature}/mash_effectsize.txt"
    return pd.read_csv(fn, sep='\t', index_col=0)\
             .rename(columns={"caudate": "Caudate",
                              "dentateGyrus": "Dentate.Gyrus",
                              "dlpfc": "DLPFC",
                              "hippocampus": "Hippocampus"})


@lru_cache()
def get_mash_obs(feature):
    fn = here("differential_analysis/tissue_comparison/_m/",feature,
              "posterior_mean_feature_4tissues.txt.gz")
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_deg(tissue, feature):
    fn = here("differential_analysis/tissue_comparison/summary_table/",
              "_m/BrainSeq_ancestry_4features_4regions.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == map_feature(feature)) &
              (df["Tissue"] == tissue.replace(".", " "))].copy()


@lru_cache()
def get_eGene(tissue, feature, lfsr):
    fn = here("eqtl_analysis/tissue_comparison/feature_summary/_m/",
              "BrainSeq_ancestry_4features_4regions.txt.gz")
    egene_df = pd.read_csv(fn, sep='\t')
    return egene_df[(egene_df["Feature"] == map_feature(feature)) &
                    (egene_df["Tissue"] == tissue.replace(".", " ")) &
                    (egene_df["lfsr"] <= lfsr)].copy()


@lru_cache()
def merge_data(feature):
    return pd.merge(get_mash_hla(feature), get_mash_obs(feature),
                    left_index=True, right_index=True,
                    suffixes=('_hla', '_obs'))


@lru_cache()
def get_sig(tissue, feature, lfsr = 0.05):
    shared_features = set(get_deg(tissue, feature).Effect) & \
        set(merge_data(feature).index)
    return merge_data(feature).loc[list(shared_features)]


@lru_cache()
def get_sig_eGene(tissue, feature, lfsr = 0.05):
    shared_features = set(get_deg(tissue, feature).Effect) & \
        set(merge_data(feature).index)
    shared_egenes = shared_features & set(get_eGene(tissue, feature, lfsr).gene_id)
    return merge_data(feature).loc[list(shared_egenes)]


@lru_cache()
def get_sig_select(tissue, feature, lfsr = 0.05):
    shared_features = set(get_deg(tissue, feature).Effect) & \
        set(merge_data(feature).index)
    shared_egenes = shared_features - set(get_eGene(tissue, feature, lfsr).gene_id)
    return merge_data(feature).loc[list(shared_egenes)]


def corr_beta(fnc, tissue, feature, lfsr = 0.05):
    return spearmanr(fnc(tissue, feature, lfsr)[f"{tissue}_hla"],
                     fnc(tissue, feature, lfsr)[f"{tissue}_obs"])


def plotNsave_corr(fnc, tissue, feature, label):
    pandas2ri.activate()
    globalenv['df'] = fnc(tissue, feature)
    globalenv['tissue'] = tissue
    globalenv['feature'] = feature
    globalenv['label'] = label
    r('''
    library(ggpubr)
    save_plot <- function(p, fn, w, h){
        for(ext in c('.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = "Observed\n(Effect Size)"
    ylab = "HLA correction\n(Effect Size)"
    fn = paste(feature, "effectsize_scatter_hla",
               tolower(tissue), label, sep=".")
    pp = ggscatter(df, x=paste0(tissue,"_obs"), y=paste0(tissue,"_hla"),
                   add="reg.line", size=1, xlab=xlab, ylab=ylab, alpha=0.5,
                   panel.labs.font=list(face="bold"),
                   add.params=list(color="blue", fill="lightgray", alpha=0.75),
                   conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                   cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                   ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=18)
    save_plot(pp, fn, 6, 6)
    ''')


def calculate_rho(fnc, feature, label):
    with open(f"rho_statistics_{feature}.{label}.log", "w") as f:
        for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
            ## Correlated effect sizes
            n = fnc(tissue, feature, 0.05).shape[0]
            rho, pval = corr_beta(fnc, tissue, feature)
            print(f"{tissue}: Predicted vs Observed (n={n}):\t"+\
                  f"r2 = {rho**2:.3f}, rho = {rho:.3f}, "+\
                  f"p-value = {pval:.3e}", file=f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Plotting correlation
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        plotNsave_corr(get_sig, tissue, args.feature, "DE")
        plotNsave_corr(get_sig_eGene, tissue, args.feature, "eGene")
        plotNsave_corr(get_sig_select, tissue, args.feature, "non_eGene")
    ## Calculate rho, DE overlap
    calculate_rho(get_sig, args.feature, "DE")
    calculate_rho(get_sig_eGene, args.feature, "eGene")
    calculate_rho(get_sig_select, args.feature, "non_eGene")
    ## Reproducibility information
    session_info.show()


if __name__ == "__main__":
    main()
