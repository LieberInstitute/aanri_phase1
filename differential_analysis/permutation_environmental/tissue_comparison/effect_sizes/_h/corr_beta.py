"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import spearmanr
from rpy2.robjects.packages import importr
from rpy2.robjects import r, globalenv, pandas2ri

@lru_cache()
def get_perm_mash(feature):
    fn = "../../summary_table/_m/"+\
        "BrainSeq_ancestry_binary_4features_4regions_allFeatures.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["feature_type"] == feature)].set_index("feature_id")


@lru_cache()
def get_perm_sig(feature):
    return get_perm_mash(feature)[(get_perm_mash(feature)["lfsr"] < 0.05)].copy()


@lru_cache()
def prepare_perm_data(feature, tissue1, tissue2):
    df1 = get_perm_sig(feature)[(get_perm_sig(feature)["region"] == tissue1)].copy()
    df2 = get_perm_sig(feature)[(get_perm_sig(feature)["region"] == tissue2)].copy()
    sig_genes = set(df1.index) | set(df2.index)
    df = get_perm_mash(feature).loc[sig_genes, ["region", "posterior_mean"]]\
                               .reset_index()
    return df[(df["region"].isin([tissue1, tissue2]))]\
        .pivot(index="feature_id", columns="region",
               values="posterior_mean")


@lru_cache()
def get_mash(feature):
    fn = here("differential_analysis/tissue_comparison/summary_table/_m",
              "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == feature)].set_index("Effect")


@lru_cache()
def get_sig(feature):
    return get_mash(feature)[(get_mash(feature)["lfsr"] < 0.05)].copy()


@lru_cache()
def prepare_data(feature, tissue1, tissue2):
    df1 = get_sig(feature)[(get_sig(feature)["Tissue"] == tissue1)].copy()
    df2 = get_sig(feature)[(get_sig(feature)["Tissue"] == tissue2)].copy()
    sig_genes = set(df1.index) | set(df2.index)
    df = get_mash(feature).loc[sig_genes, ["Tissue", "posterior_mean"]]\
                          .reset_index()
    return df[(df["Tissue"].isin([tissue1, tissue2]))]\
        .pivot(index="Effect", columns="Tissue",
               values="posterior_mean")


@lru_cache()
def prepare_environment(feature, tissue):
    df1 = get_mash(feature).loc[(get_mash(feature)["Tissue"] == tissue),
                                ["lfsr", "posterior_mean"]]\
                           .rename(columns={"posterior_mean":"continous"}).copy()
    df2 = get_perm_mash(feature).loc[(get_perm_mash(feature)["region"] == tissue),
                                     ["lfsr", "posterior_mean"]]\
                                .rename(columns={"posterior_mean":"binary"}).copy()
    sig_genes = set(df1[(df1["lfsr"] < 0.05)].index) | \
        set(df2[(df2["lfsr"] < 0.05)].index)
    df = pd.merge(df1, df2, left_index=True, right_index=True)\
           .drop(["lfsr_x", "lfsr_y"], axis=1)\
           .reset_index()
    return df[(df["index"].isin(sig_genes))].copy()


@lru_cache()
def prepare_environment_shared(feature, tissue):
    df1 = get_mash(feature).loc[(get_mash(feature)["Tissue"] == tissue),
                                ["lfsr", "posterior_mean"]]\
                           .rename(columns={"posterior_mean":"continous"}).copy()
    df2 = get_perm_mash(feature).loc[(get_perm_mash(feature)["region"] == tissue),
                                     ["lfsr", "posterior_mean"]]\
                                .rename(columns={"posterior_mean":"binary"}).copy()
    sig_genes = set(df1[(df1["lfsr"] < 0.05)].index) & \
        set(df2[(df2["lfsr"] < 0.05)].index)
    df = pd.merge(df1, df2, left_index=True, right_index=True)\
           .drop(["lfsr_x", "lfsr_y"], axis=1)\
           .reset_index()
    return df[(df["index"].isin(sig_genes))].copy()


def corr_beta(tissue1, tissue2):
    df = prepare_perm_data("Gene", tissue1, tissue2)
    return spearmanr(df[tissue1], df[tissue2])


def corr_beta_environment(fnc, tissue, feature):
    df = fnc(feature, tissue)
    return spearmanr(df["continous"], df["binary"])


def plotNsave_corr(tissue1, tissue2):
    pandas2ri.activate()
    globalenv['df'] = prepare_perm_data("Gene", tissue1, tissue2)
    globalenv['tissue1'] = tissue1
    globalenv['tissue2'] = tissue2
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = paste0("Effect Size\n(", tissue1, ")")
    ylab = paste0("Effect Size\n(", tissue2, ")")
    fn = paste("effectsize_scatter", gsub(" ", "_", tissue1),
               gsub(" ", "_", tissue2), sep="_")
    pp = ggpubr::ggscatter(df, x=tissue1, y=tissue2, add="reg.line", size=1,
                           xlab=xlab, ylab=ylab,panel.labs.font=list(face="bold"),
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15), ncol=4)
    save_plot(pp, fn, 6, 6)
    ''')


def plotNsave_corr_environment(fnc, tissue, label, feature):
    pandas2ri.activate()
    globalenv['label'] = f"{label}_{feature.lower()}"
    globalenv['df'] = fnc(feature, tissue)
    globalenv['tissue'] = tissue
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = paste0("Effect Size\n(Black Americans)")
    ylab = paste0("Effect Size\n(Black & White Americans)")
    fn = paste(label, "environmental_effectsize_scatter",
               gsub(" ", "_", tissue), sep="_")
    pp = ggpubr::ggscatter(df, x="continous", y="binary", add="reg.line", size=1,
                           alpha=0.3, xlab=xlab, ylab=ylab,
                           panel.labs.font=list(face="bold"),
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15), ncol=4)
    save_plot(pp, fn, 6, 6)
    ''')


def main():
    ## Comparison
    config = {'either': prepare_environment,
              'shared': prepare_environment_shared}
    with open("rho_statistics_environment.log", "w") as f:
        for feature in ["Gene", "Transcript", "Exon", "Junction"]:
            print(feature.upper(), file=f)
            for label in ["either", "shared"]:
                print(label.upper(), file=f)
                for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
                    ## Generate figure
                    plotNsave_corr_environment(config[label], tissue, label, feature)
                    ## Correlated effect sizes
                    rho, pval = corr_beta_environment(config[label], tissue, feature)
                    print(f"Environmental ({tissue}):\t\t r2 ~ {rho**2:.3f}\t"+\
                          f"rho > {rho:.3f}\tp-value < {pval:.1e}", file=f)
            
    ## Calculate rho
    with open("rho_statistics.log", "w") as f:
        for tissue1 in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            for tissue2 in ["Dentate Gyrus", "DLPFC", "Hippocampus"]:
                if tissue1 != tissue2:
                    ## Generate figure
                    plotNsave_corr(tissue1, tissue2)
                    ## Correlated effect sizes
                    rho, pval = corr_beta(tissue1, tissue2)
                    print("%s vs %s:\t\t\t rho > %.3f, p-value < %.1e" %
                          (tissue1, tissue2, rho, pval), file=f)


if __name__ == "__main__":
    main()
