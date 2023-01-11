"""
This script preforms the Fisher's pairwise enrichment analysis
between brain regions.
"""
import numpy as np
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

@lru_cache()
def get_public_deg():
    fn = here("input/public/nedelec_immune_bulk.xlsx")
    return pd.read_excel(fn, sheet_name="B) pop-DE", skiprows=3)


@lru_cache()
def get_mash_lfsr(tissue):
    # get lfsr
    return pd.read_csv("../../_m/genes/lfsr_feature_4tissues.txt.gz", sep='\t')\
             .loc[:, ["Effect", tissue]]\
             .rename(columns={tissue: "lfsr"})


@lru_cache()
def get_mash_es(tissue):
    # get effect size
    return pd.read_csv("../../_m/genes/posterior_mean_feature_4tissues.txt.gz", sep='\t')\
             .loc[:, ["Effect", tissue]]\
             .rename(columns={tissue: "posterior_mean"})


@lru_cache()
def annotate_mash(tissue):
    df = get_mash_lfsr(tissue)\
        .merge(get_mash_es(tissue), on="Effect")
    df["ensembl_gene_id"] = df.Effect.str.replace("\\..*", "", regex=True)
    return df


@lru_cache()
def merge_data(tissue):
    """
    Merge data to get shared features to check for enrichment.
    """
    return annotate_mash(tissue).merge(get_public_deg(), on="ensembl_gene_id")


def cal_fishers(tissue, lab, direction):
    df = merge_data(tissue)
    if direction == 'Up':
        df = df[(df['posterior_mean'] > 0)].copy()
    elif direction == 'Down':
        df = df[(df['posterior_mean'] < 0)].copy()
    else:
        df = df
    table = [[sum((df["lfsr"]<0.05) & ((df[lab]<0.05))),
              sum((df["lfsr"]<0.05) & ((df[lab]>=0.05)))],
             [sum((df["lfsr"]>=0.05) & ((df[lab]<0.05))),
              sum((df["lfsr"]>=0.05) & ((df[lab]>=0.05)))]]
    # print(table)
    return fisher_exact(table)


def enrichment_loop():
    comp_list = {"NI_fdr": "Non-infected", "L_fdr": "Listeria",
                 "S_fdr": "Salmonella"}
    or_lt = []; pval_lt = []; tissue_lt = []; comp_lt = []; dir_lt = []; fdr_lt = []
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        pvals = []
        for direction in ["Up", "Down", "All"]:
            for comp in ["NI_fdr", "L_fdr", "S_fdr"]:
                ## Correlated effect sizes
                odds, pval = cal_fishers(tissue, comp, direction)
                or_lt.append(odds); pvals.append(pval); tissue_lt.append(tissue);
                comp_lt.append(comp_list[comp]); dir_lt.append(direction)
        fdr = multipletests(pvals, method='fdr_bh')[1]
        pval_lt = np.concatenate((pval_lt, pvals))
        fdr_lt = np.concatenate((fdr_lt, fdr))
    return pd.DataFrame({"Tissue": tissue_lt, "Conditions": comp_lt, "OR": or_lt,
                         "P-value": pval_lt, "FDR": fdr_lt, "Direction": dir_lt})


def main():
    ## Enrichment loop and save results
    enrichment_loop()\
        .to_csv("nedelec_immune_enrichment_analysis_DEGs.tsv",
                sep='\t', index=False)

    ## Calculate rho
    with open("fishers_enrichment.log", "w") as f:
        for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
            for label in ["NI_fdr", "L_fdr", "S_fdr"]:
                for direction in ["Up", "Down", "All"]:
                    ## Correlated effect sizes
                    odds, pval = cal_fishers(tissue, label, direction)
                    print("%s:\t%s vs %s:\tOdds Ratio > %.3f, p-value < %.1e" %
                          (direction, tissue, label, odds, pval), file=f)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
