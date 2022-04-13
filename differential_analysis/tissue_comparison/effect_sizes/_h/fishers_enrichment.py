"""
This script preforms the Fisher's pairwise enrichment analysis
between brain regions.
"""
import pandas as pd
from functools import lru_cache
from scipy.stats import fisher_exact

@lru_cache()
def get_mash():
    # get effect size
    fn = "../../_m/genes/lfsr_feature_4tissues.txt.gz"
    return pd.read_csv(fn, sep='\t', index_col=0)


def cal_fishers(tissue1, tissue2):
    table = [[sum((get_mash()[tissue1]<0.05) & ((get_mash()[tissue2]<0.05))),
              sum((get_mash()[tissue1]<0.05) & ((get_mash()[tissue2]>=0.05)))],
             [sum((get_mash()[tissue1]>=0.05) & ((get_mash()[tissue2]<0.05))),
              sum((get_mash()[tissue1]>=0.05) & ((get_mash()[tissue2]>=0.05)))]]
    print(table)
    return fisher_exact(table, alternative='greater')


def main():
    ## Calculate rho
    with open("fishers_enrichment.log", "w") as f:
        for tissue1 in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
            for tissue2 in ["Dentate.Gyrus", "DLPFC", "Hippocampus"]:
                if tissue1 != tissue2:
                    ## Correlated effect sizes
                    odds, pval = cal_fishers(tissue1, tissue2)
                    print("%s vs %s:\t\t\t  Odds Ratio > %.3f, p-value < %.1e" %
                          (tissue1, tissue2, odds, pval), file=f)


if __name__ == "__main__":
    main()
