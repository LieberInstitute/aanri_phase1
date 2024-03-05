"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import argparse
import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import spearmanr

def map_feature(feature):
    return {"genes": "Gene", "transcripts": "Transcript",
            "exons": "Exon", "junctions": "Junction"}[feature]


@lru_cache()
def get_mash_pred(feature):
    # get effect size
    fn = f"../../_m/{feature}/posterior_mean_feature_4tissues.txt.gz"
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_mash_obs(feature):
    fn = here(f"differential_analysis/tissue_comparison/_m/{feature}/",
              "posterior_mean_feature_4tissues.txt.gz")
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_deg(tissue, feature):
    fn = here("differential_analysis/tissue_comparison/summary_table/",
              "_m/BrainSeq_ancestry_4features_4regions.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == map_feature(feature)) &
              (df["Tissue"] == tissue.replace(".", " "))].copy()


@lru_cache()
def get_eGene(tissue, feature, lfsr):
    fn = here("eqtl_analysis/tissue_comparison/feature_summary/_m/",
              "BrainSeq_ancestry_4features_4regions.txt.gz")
    egene_df = pd.read_csv(fn, sep='\t')
    return egene_df[(egene_df["Feature"] == map_feature(feature)) &
                    (egene_df["Tissue"] == tissue.replace(".", " ")) &
                    (egene_df["lfsr"] <= lfsr)].copy()


@lru_cache()
def merge_data(feature):
    return pd.merge(get_mash_pred(feature), get_mash_obs(feature),
                    left_index=True, right_index=True,
                    suffixes=('_pred', '_obs'))


@lru_cache()
def get_sig(tissue, feature, lfsr = 0.05):
    shared_features = set(get_deg(tissue, feature).Effect) & \
        set(merge_data(feature).index)
    return merge_data(feature).loc[list(shared_features)]


@lru_cache()
def get_sig_eGene(tissue, feature, lfsr = 0.05):
    shared_features = set(get_deg(tissue, feature).Effect) & \
        set(merge_data(feature).index)
    shared_egenes = shared_features & set(get_eGene(tissue, feature, lfsr).gene_id)
    return merge_data(feature).loc[list(shared_egenes)]


@lru_cache()
def get_sig_select(tissue, feature, lfsr = 0.05):
    shared_features = set(get_deg(tissue, feature).Effect) & \
        set(merge_data(feature).index)
    shared_egenes = shared_features - set(get_eGene(tissue, feature, lfsr).gene_id)
    return merge_data(feature).loc[list(shared_egenes)]


def corr_beta(fnc, tissue, feature, lfsr = 0.05):
    return spearmanr(fnc(tissue, feature, lfsr)[f"{tissue}_pred"],
                     fnc(tissue, feature, lfsr)[f"{tissue}_obs"])


def calculate_rho(fnc, feature, label):
    with open(f"rho_statistics_{feature}.{label}.log", "w") as f:
        for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
            ## Correlated effect sizes
            n = fnc(tissue, feature, 0.05).shape[0]
            rho, pval = corr_beta(fnc, tissue, feature)
            print(f"{tissue}: Predicted vs Observed (n={n}):\t"+\
                  f"r2 = {rho**2:.3f}, rho = {rho:.3f}, "+\
                  f"p-value = {pval:.3e}", file=f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Calculate rho, DE overlap
    calculate_rho(get_sig, args.feature, "DE")
    calculate_rho(get_sig_eGene, args.feature, "eGene")
    calculate_rho(get_sig_select, args.feature, "non_eGene")
    ## Reproducibility information
    session_info.show()


if __name__ == "__main__":
    main()
