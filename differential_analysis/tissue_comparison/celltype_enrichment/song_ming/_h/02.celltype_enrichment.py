## Cell type enrichment analysis

import numpy as np
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

@lru_cache()
def get_celltype_genes(dataset):
    fn = f'{dataset.lower()}.single_cell.tsv'
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_ancestry_degs():
    fn = "../../../summary_table/_m/"+\
        "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == "Gene")].copy()


@lru_cache()
def subset_region(region):
    df = get_ancestry_degs()
    return df.loc[(df["Tissue"] == region),
                  ["gencodeID", "Symbol",
                   "lfsr", "posterior_mean"]]\
             .dropna(subset=["Symbol"])


@lru_cache()
def merge_dataframe(celltype, tissue, dataset):
    return subset_region(tissue)\
        .merge(get_celltype_genes(dataset).loc[:, [celltype]],
               left_on='Symbol', right_index=True)


def cal_fishers_direction(celltype, direction, tissue, dataset):
    df = merge_dataframe(celltype, tissue, dataset)
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


def calculate_enrichment(dataset):
    region_lt = []; dir_lt = []; fdr_lt = [];
    ct_lt = []; pval_lt = []; oddratio_lt = []
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        pvals = []
        for direction in ['Up', 'Down', 'All']:
            for celltype in get_celltype_genes(dataset).columns.unique():
                odd_ratio, pval = cal_fishers_direction(celltype, direction, tissue, dataset)
                ct_lt.append(celltype); pvals.append(pval); region_lt.append(tissue)
                oddratio_lt.append(odd_ratio); dir_lt.append(direction)
        _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
        pval_lt = np.concatenate((pval_lt, pvals))
        fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, 'Cell_type': ct_lt, 'OR': oddratio_lt,
                         'PValue': pval_lt, "FDR": fdr_lt, 'Direction': dir_lt})


def main():
    ## Calculate enrichment
    for dataset in ["microglia", "astrocyte", "oligo"]:
        print(dataset)
        calculate_enrichment(dataset)\
            .to_csv(f'celltype_enrichment_analysis.{dataset}.txt',
                    sep='\t', index=False)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
    
