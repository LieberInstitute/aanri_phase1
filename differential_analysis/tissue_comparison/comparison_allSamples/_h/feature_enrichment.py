# Examine enrichment in psychiatric disorders TWAS and DEGs
import numpy as np
import pandas as pd
from os import environ
from functools import lru_cache
from scipy.stats import spearmanr
from scipy.stats import fisher_exact
from rpy2.robjects import r, globalenv, pandas2ri
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '32'

@lru_cache()
def get_background(tissue, feature):
    f_dict = {"Gene":"genes", "Transcript": "transcripts",
              "Exon": "exons", "Junction": "junctions"}
    fn = "../../../../internal_replication/tissue_comparison/"+\
        "mash/_m/%s/lfsr_feature_4tissues.txt.gz" % f_dict[feature]
    return pd.read_csv(fn, sep='\t')\
             .loc[:, ["Effect", tissue.replace(" ", ".")]]


@lru_cache()
def get_bg_aa(tissue, feature):
    f_dict = {"Gene":"genes", "Transcript": "transcripts",
              "Exon": "exons", "Junction": "junctions"}
    fn = "../../_m/%s/lfsr_feature_4tissues.txt.gz" % f_dict[feature]
    return pd.read_csv(fn, sep='\t')\
             .loc[:, ["Effect", tissue.replace(" ", ".")]]


@lru_cache()
def get_mash(tissue, feature):
    fn = "../../../../internal_replication/tissue_comparison/mash/"+\
        "summary_table/_m/BrainSeq_ancestry_4features_4regions.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Tissue"] == tissue) & (df["Type"] == feature)].copy()


@lru_cache()
def get_mash_aa(tissue, feature):
    fn = "../../summary_table/_m/BrainSeq_ancestry_4features_4regions.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Tissue"] == tissue) & (df["Type"] == feature)].copy()


@lru_cache()
def get_es_all(tissue, feature):
    f_dict = {"Gene":"genes", "Transcript": "transcripts",
              "Exon": "exons", "Junction": "junctions"}
    fn = "../../../../internal_replication/tissue_comparison/mash/_m/"+\
        "%s/posterior_mean_feature_4tissues.txt.gz" % f_dict[feature]
    return pd.read_csv(fn, sep='\t')\
             .loc[:, ["Effect", tissue.replace(" ", ".")]]


@lru_cache()
def prepare_data(tissue, feature):
    #feature = "Junction"
    return get_mash_aa(tissue, feature)\
        .merge(get_es_all(tissue, feature)\
               .rename(columns={tissue.replace(" ", "."): "ALL"}),
               on="Effect")


def fet(tissue, feature):
    """
    Calculates Fisher's Exact test (fet) with sets a and b in universe u.
    Inputs are sets.
    """
    u_a = set(get_background(tissue, feature).Effect)
    u_b = set(get_bg_aa(tissue, feature).Effect)
    a = set(get_mash(tissue, feature).Effect)
    b = set(get_mash_aa(tissue, feature).Effect)
    yes_a = u_a.intersection(a)
    yes_b = u_b.intersection(b)
    no_a = u_a - a
    no_b = u_b - b
    m = [[len(yes_a.intersection(yes_b)), len(no_a.intersection(yes_b))],
         [len(yes_a.intersection(no_b)), len(no_a.intersection(no_b))]]
    #print(m)
    return fisher_exact(m)


def enrichment_loop(feature):
    or_lt = []; pval_lt = []; tissue_lt = []
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        oddratio, pvals = fet(tissue, feature)
        or_lt.append(oddratio); pval_lt.append(pvals);
        tissue_lt.append(tissue);
    fdr = multipletests(pval_lt, method='fdr_bh')[1]
    return pd.DataFrame({"Tissue": tissue_lt, "Feature": feature,
                         "OR": or_lt, "P-value": pval_lt, "FDR": fdr})


def corr_beta(tissue, feature):
    return spearmanr(prepare_data(tissue, feature)["posterior_mean"],
                     prepare_data(tissue, feature)["ALL"])


def plotNsave_corr(tissue, feature):
    pandas2ri.activate()
    globalenv['df'] = prepare_data(tissue, feature).loc[:, ["posterior_mean", "ALL"]]
    globalenv['tissue'] = tissue
    globalenv['feature'] = feature
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.png', '.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    ylab = paste0("Effect Size\n(All Samples)")
    xlab = paste0("Effect Size\n(AA only)")
    fn = paste("effectsize_scatter", gsub(" ", "_", tissue), feature, sep="_")
    pp = ggpubr::ggscatter(df, x="posterior_mean", y="ALL", add="reg.line", size=1,
                           xlab=xlab, ylab=ylab,panel.labs.font=list(face="bold"),
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15), ncol=4)
    save_plot(pp, fn, 6, 6)
    ''')


def main():
    ## Enrichment
    df = pd.DataFrame()
    for feature in ["Gene", "Transcript", "Exon", "Junction"]:
        df = pd.concat([df, enrichment_loop(feature)], axis=0)
    ## Save enrichment
    df.to_csv("enrichment_analysis_all_vs_AAonly.tsv",
              sep='\t', index=False)
    with open("rho_statistics.log", "w") as f:
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            print("%s:\n" % tissue, file=f)
            for feature in ["Gene", "Transcript", "Exon", "Junction"]:
                ## Generate figure
                plotNsave_corr(tissue, feature)
                ## Correlated effect sizes
                rho, pval = corr_beta(tissue, feature)
                print("%s:\t rho > %.3f, p-value < %.1e" %
                      (feature, rho, pval), file=f)


if __name__ == '__main__':
    main()
