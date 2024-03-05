# Examine enrichment of elastic net (r2 > 0.01)

import numpy as np
import pandas as pd
import session_info
from os import environ
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '5'

def map_feature(feature):
    return {"genes":"Gene", "transcripts":"Transcript",
            "exons":"Exon", "junctions":"Junction"}[feature]


@lru_cache()
def get_metrics():
    return pd.read_csv("../_m/prediction_metrics_summary.txt.gz",
                       sep='\t', compression="gzip")


@lru_cache()
def get_background(tissue, feature):
    new_tissue = tissue.replace(" ", "_")
    return get_metrics()[(get_metrics()["region"] == new_tissue) &
                         (get_metrics()["type"] == map_feature(feature))]


@lru_cache()
def get_predicted(tissue, feature, r2=0.01):
    new_tissue = tissue.replace(" ", "_")
    return get_metrics()[(get_metrics()["region"] == new_tissue) &
                         (get_metrics()["type"] == map_feature(feature)) &
                         (get_metrics()["test_score_r2_median"] > r2)].copy()


@lru_cache()
def get_de(tissue, feature):
    fn = here("differential_analysis/tissue_comparison/summary_table/",
              "_m/BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Tissue"] == tissue) &
              (df["Type"] == map_feature(feature))].copy()


@lru_cache()
def get_ancestry_deg(tissue, feature, direction):
    df = get_de(tissue, feature)
    if direction=="all":
        return df[(df["lfsr"] < 0.05)].copy()
    elif direction=="up":
        return df[(df["lfsr"] < 0.05) & (df["posterior_mean"] > 0)].copy()
    else:
        return df[(df["lfsr"] < 0.05) & (df["posterior_mean"] < 0)].copy()


def cal_enrichment(tissue, feature, direction, r2=0.01):
    """
    Calculates Fisher's Exact test.
    Inputs: brainr region and DE direction of effect.
    """
    universe = set(get_background(tissue, feature).feature) | \
        set(get_de(tissue, feature).Effect)
    de_set   = set(get_ancestry_deg(tissue, feature, direction).Effect)
    pred_set = set(get_predicted(tissue, feature, r2).feature)
    yes_de   = universe.intersection(de_set)
    yes_pred = universe.intersection(pred_set)
    no_de    = universe - de_set
    no_pred  = universe - pred_set
    m = [[len(yes_de.intersection(yes_pred)),
          len(no_de.intersection(yes_pred))],
         [len(yes_de.intersection(no_pred)),
          len(no_de.intersection(no_pred))]]
    odds, pval = fisher_exact(m)
    n = len(yes_de.intersection(yes_pred))
    frac_overlap = n / len(yes_pred)
    return odds, pval, frac_overlap, n


def enrichment_loop():
    dir_dict = {"all": "All", "up": "Decreased in AA",
                "down": "Increased in AA"}
    dt = pd.DataFrame()
    for direction in ["all", "up", "down"]:
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            for feature in ["genes", "transcripts", "exons", "junctions"]:
                or_lt = []; pval_lt = []; tissue_lt = []; r2_lt = [];
                frac_lt = []; n_lt = [ ];
                for r2 in np.append(0.01, np.arange(0.05, 0.55, 0.05)):
                    oddratio, pvals, frac, n = cal_enrichment(tissue, feature,
                                                              direction, r2)
                    or_lt.append(oddratio); pval_lt.append(pvals);
                    tissue_lt.append(tissue); r2_lt.append(r2);
                    frac_lt.append(frac); n_lt.append(n)
                fdr = multipletests(pval_lt, method='fdr_bh')[1]
                dtx = pd.DataFrame({"Tissue": tissue_lt,
                                    "Feature": map_feature(feature),
                                    "Direction": dir_dict[direction],
                                    "R2_threshold": r2_lt, "OR": or_lt,
                                    "P-value": pval_lt, "FDR": fdr,
                                    "Fraction_Overlap": frac_lt,
                                    "N_Overlap": n_lt})
                dt = pd.concat([dt, dtx], axis=0)
    return dt


def main():
    ## Enrichment analysis
    dt = enrichment_loop()
    ## Save enrichment
    dt.to_csv("predicted_r2_enrichment_fishers.tsv", sep='\t', index=False)
    session_info.show()


if __name__ == '__main__':
    main()
