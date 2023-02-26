"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

@lru_cache()
def get_local_mash(feature):
    fn = "../../summary_table/_m/"+\
        "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["feature_type"] == feature)].set_index("feature_id")


@lru_cache()
def get_mash(feature):
    fn = here("differential_analysis/tissue_comparison/summary_table/_m",
              "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == feature)].set_index("Effect")


@lru_cache()
def prepare_data(feature, tissue):
    df1 = get_mash(feature).loc[(get_mash(feature)["Tissue"] == tissue),
                                ["lfsr", "posterior_mean"]]\
                           .rename(columns={"posterior_mean":"continous"}).copy()
    df2 = get_local_mash(feature).loc[(get_local_mash(feature)["region"] == tissue),
                                     ["lfsr", "posterior_mean"]]\
                                .rename(columns={"posterior_mean":"local"}).copy()
    return pd.merge(df1, df2, left_index=True, right_index=True,
                    suffixes=["_global", "_local"])


def cal_fishers_feature(feature, tissue, direction):
    df = prepare_data(feature, tissue)
    if direction == 'Up': ## using global ancestry
        df = df[(df['continous'] > 0)].copy()
    elif direction == 'Down':
        df = df[(df['continous'] < 0)].copy()
    else:
        df = df
    table = [[np.sum((df['lfsr_global']<=0.05) & (df['lfsr_local']<=0.05)),
              np.sum((df['lfsr_global']<=0.05) & (df['lfsr_local']>0.05))],
             [np.sum((df['lfsr_global']>0.05)  & (df['lfsr_local']<=0.05)),
              np.sum((df['lfsr_global']>0.05)  & (df['lfsr_local']>0.05))]]
    print(f"Overlap {direction}: {table[0][0]} ({table[0][0] / np.sum(table[0]):.1%})")
    return fisher_exact(table)


def calculate_enrichment():
    region_lt = []; feat_lt = []; fdr_lt = []; pval_lt = [];
    dir_lt = []; oddratio_lt = []
    for feature in ['Gene', 'Transcript', 'Exon', 'Junction']:
        print(feature)
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            print(tissue)
            pvals = []
            for direction in ['Up', 'Down', 'All']:
                odd_ratio, pval = cal_fishers_feature(feature, tissue, direction)
                pvals.append(pval); region_lt.append(tissue)
                oddratio_lt.append(odd_ratio); feat_lt.append(feature)
                dir_lt.append(direction)
            _, fdr = fdrcorrection(pvals) # FDR correction per feature
            pval_lt = np.concatenate((pval_lt, pvals))
            fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate DataFrame
    return pd.DataFrame({'Region': region_lt, 'Feature': feat_lt,
                         'OR': oddratio_lt, 'PValue': pval_lt, "FDR": fdr_lt,
                         'Direction': dir_lt})


def main():
    ## Global vs Local ancestry enrichment / overlap
    calculate_enrichment()\
        .to_csv('enrichment_analysis.txt', sep='\t', index=False)

    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
