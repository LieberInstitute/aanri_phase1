## Cell type enrichment analysis

import numpy as np
import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

@lru_cache()
def get_mhc_genes():
    fn = '../_h/xmhc_genes_classification.csv'
    return pd.read_csv(fn, sep=',', usecols=[3,5,6])


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
    # df["subregion"] = df["subregion"].fillna("NA")
    # df["gene_cluster"] = df["gene_cluster"].fillna("NA")
    return df[(df["Tissue"] == tissue)].copy()


def cal_fishers_cluster(direction, tissue, gcluster):
    df = merge_dataframe(tissue)
    if direction == 'Up':
        df = df[(df['posterior_mean'] > 0)].copy()
    elif direction == 'Down':
        df = df[(df['posterior_mean'] < 0)].copy()
    else:
        df = df
    table = [[np.sum((df['lfsr']<0.05) & (df["gene_cluster"] == gcluster)),
              np.sum((df['lfsr']<0.05) & (df["gene_cluster"] != gcluster))],
             [np.sum((df['lfsr']>0.05) & (df["gene_cluster"] == gcluster)),
              np.sum((df['lfsr']>0.05) & (df["gene_cluster"] != gcluster))]]
    return fisher_exact(table)


def cal_fishers_subregion(direction, tissue, subregion):
    df = merge_dataframe(tissue)
    if direction == 'Up':
        df = df[(df['posterior_mean'] > 0)].copy()
    elif direction == 'Down':
        df = df[(df['posterior_mean'] < 0)].copy()
    else:
        df = df
    table = [[np.sum((df['lfsr']<0.05) & (df["subregion"] == subregion)),
              np.sum((df['lfsr']<0.05) & (df["subregion"] != subregion))],
             [np.sum((df['lfsr']>0.05) & (df["subregion"] == subregion)),
              np.sum((df['lfsr']>0.05) & (df["subregion"] != subregion))]]
    return fisher_exact(table)


def calculate_enrichment(label, lnames, fnc):
    region_lt = []; dir_lt = []; fdr_lt = []; pval_lt = []; oddratio_lt = []
    label_lt = [];
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        pvals = []
        for xlabel in lnames:
            for direction in ['Up', 'Down', 'All']:
                odd_ratio, pval = fnc(direction, tissue, xlabel)
                pvals.append(pval); region_lt.append(tissue)
                oddratio_lt.append(odd_ratio); dir_lt.append(direction)
                label_lt.append(xlabel)
        _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
        pval_lt = np.concatenate((pval_lt, pvals))
        fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Region': region_lt, 'Direction': dir_lt,
                         label: label_lt, 'OR': oddratio_lt,
                         'PValue': pval_lt, "FDR": fdr_lt})


def main():
    # Gene clusters
    clusters = get_mhc_genes().dropna(subset="gene_cluster")\
                              .gene_cluster.unique()
    dx = calculate_enrichment("Gene_Cluster", clusters,
                              cal_fishers_cluster)
    dx.to_csv('xmhc_enrichment_analysis.gene_cluster.txt',
              sep='\t', index=False)
    # Subregion
    subregions = get_mhc_genes().dropna(subset="subregion")\
                                .subregion.unique()
    dx = calculate_enrichment("Subregion", subregions,
                              cal_fishers_subregion)
    dx.to_csv('xmhc_enrichment_analysis.subregion.txt',
              sep='\t', index=False)
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
    
