## Cell type enrichment analysis

import numpy as np
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

@lru_cache()
def get_mhc_genes():
    df = pd.read_csv('../_m/xmhc_genes.csv', sep=',', usecols=[3])
    df["mhc"] = 1
    return df


@lru_cache()
def get_degs():
    fn = "../../summary_table/_m/BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz"
    df = pd.read_csv(fn, sep="\t")
    return df[(df["Type"] == "Gene")].copy()


@lru_cache()
def merge_dataframe(tissue):
    df = pd.merge(get_mhc_genes(), get_degs(),
                  left_on='gene_id', right_on="gencodeID",
                  how="right")
    df["mhc"] = df["mhc"].fillna(0)
    return df[(df["Tissue"] == tissue)].copy()


def cal_fishers_direction(direction, tissue):
    df = merge_dataframe(tissue)
    if direction == 'Up':
        df = df[(df['posterior_mean'] > 0)].copy()
    elif direction == 'Down':
        df = df[(df['posterior_mean'] < 0)].copy()
    else:
        df = df
    table = [[np.sum((df['lfsr']<0.05) & (df["mhc"] == 1)),
              np.sum((df['lfsr']<0.05) & (df["mhc"] == 0))],
             [np.sum((df['lfsr']>0.05) & (df["mhc"] == 1)),
              np.sum((df['lfsr']>0.05) & (df["mhc"] == 0))]]
    print(table)
    return fisher_exact(table)


def calculate_enrichment():
    region_lt = []; dir_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        pvals = []
        for direction in ['Up', 'Down', 'All']:
            odd_ratio, pval = cal_fishers_direction(direction, tissue)
            pvals.append(pval); region_lt.append(tissue)
            oddratio_lt.append(odd_ratio); dir_lt.append(direction)
        _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
        pval_lt = np.concatenate((pval_lt, pvals))
        fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Region': region_lt, 'Direction': dir_lt,
                         'OR': oddratio_lt,'PValue': pval_lt,"FDR": fdr_lt})


def main():
    calculate_enrichment()\
        .to_csv('xmhc_enrichment_analysis.txt', sep='\t', index=False)


if __name__ == "__main__":
    main()
    
