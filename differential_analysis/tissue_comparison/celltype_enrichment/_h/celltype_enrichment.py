## Cell type enrichment analysis

import numpy as np
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

@lru_cache()
def get_celltype_genes():
    fn = here('input/celltypes/_m/Zeisel_single_cell.tsv')
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_annotation():
    base_loc = here("input/text_files_counts")
    config = {
        "genes": here("input/text_files_counts/_m",
                      "caudate/gene_annotation.tsv"),
        "transcripts": here("input/text_files_counts/_m",
                            "caudate/tx_annotation.tsv"),
        "exons": here("input/text_files_counts/_m",
                      "caudate/exon_annotation.tsv"),
        "junctions": here("input/text_files_counts/_m",
                          "caudate/jxn_annotation.tsv"),
    }
    return pd.read_csv(config["genes"], sep='\t')\
             .loc[:, ["names", "seqnames", "start",
                      "end", "Symbol", "gencodeID"]]


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
    return get_mash_lfsr(tissue)\
        .merge(get_mash_es(tissue), on="Effect")\
        .merge(get_annotation(), left_on="Effect", right_on="names")\
        .drop(["names"], axis=1).dropna(subset=["Symbol"])


@lru_cache()
def merge_dataframe(celltype, tissue):
    return annotate_mash(tissue).merge(get_celltype_genes().loc[:, [celltype]],
                                       left_on='Symbol', right_index=True)


def cal_fishers_direction(celltype, direction, tissue):
    df = merge_dataframe(celltype, tissue)
    if direction == 'Up':
        df = df[(df['posterior_mean'] > 0)].copy()
    elif direction == 'Down':
        df = df[(df['posterior_mean'] < 0)].copy()
    else:
        df = df
    
    table = [[np.sum((df['lfsr']<0.05) & (df[celltype] == 1)), 
              np.sum((df['lfsr']<0.05) & (df[celltype] == 0))],
             [np.sum((df['lfsr']>0.05) & (df[celltype] == 1)), 
              np.sum((df['lfsr']>0.05) & (df[celltype] == 0))]]
    print(table)
    return fisher_exact(table)


def calculate_enrichment():
    region_lt = []; dir_lt = []; fdr_lt = []; ct_lt = []; pval_lt = []; oddratio_lt = []
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        pvals = []
        for direction in ['Up', 'Down', 'All']:
            for celltype in get_celltype_genes().columns.unique():
                odd_ratio, pval = cal_fishers_direction(celltype, direction, tissue)
                ct_lt.append(celltype); pvals.append(pval); region_lt.append(tissue)
                oddratio_lt.append(odd_ratio); dir_lt.append(direction)
        _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
        pval_lt = np.concatenate((pval_lt, pvals))
        fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, 'Cell_type': ct_lt, 'OR': oddratio_lt,
                         'PValue': pval_lt, "FDR": fdr_lt, 'Direction': dir_lt})


def main():
    calculate_enrichment().to_csv('celltype_enrichment_analysis.txt', sep='\t', index=False)


if __name__ == "__main__":
    main()
    
