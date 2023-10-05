"""
This script selects genes identified in Black only compared with
permutation Black v. white analysis.
"""
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache

@lru_cache()
def get_environmental_degs(feature, region):
    fn = "../../summary_table/_m/BrainSeq_ancestry_binary_4features_4regions.txt.gz"
    df = pd.read_csv(fn, sep="\t")
    df["environ"] = 1
    return df.loc[(df["feature_type"] == feature) & (df["region"] == region),
                  ["feature_id", "gene_id", "environ"]].copy()


@lru_cache()
def get_ancestry_degs(feature, region):
    fn = here("differential_analysis/tissue_comparison/summary_table",
              "_m/BrainSeq_ancestry_4features_4regions.txt.gz")
    df = pd.read_csv(fn, sep="\t")
    df["ancestry"] = 1
    return df.loc[(df["Type"] == feature) & (df["Tissue"] == region),
                  ["Effect", "ancestry"]].copy()


@lru_cache()
def de_comparison(feature, region):
    genes = set(get_environmental_degs(feature, region).feature_id) | \
        set(get_ancestry_degs(feature, region).Effect)
    return pd.DataFrame({"feature_id":list(genes)})\
             .merge(get_environmental_degs(feature, region), on="feature_id", how="left")\
             .merge(get_ancestry_degs(feature, region),
                    left_on="feature_id", right_on="Effect", how="left")\
             .drop(["Effect"], axis=1).fillna(0)


def main():
    # Merge data
    dt = pd.DataFrame()
    for feature_type in ["Gene", "Transcript", "Exon", "Junction"]:
        for region in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            df = de_comparison(feature_type, region)
            df["region"] = region; df["feature_type"] = feature_type
            dt = pd.concat([dt, df], axis=0)
    dt.to_csv("BrainSeq_ancestry_DE_comparison.txt.gz", sep='\t', index=False)
    # Summarize results: ancestry specific (Black only)
    ancestry_specific = dt[(dt["ancestry"] == 1) & (dt["environ"] != 1)].copy()
    print("Ancestry specific:")
    print(ancestry_specific.groupby(["region", "feature_type"]).size())
    # Summary: potential environmental DE
    environ_specific = dt[(dt["ancestry"] != 1) & (dt["environ"] == 1)].copy()
    print("Environmental specific:")
    print(environ_specific.groupby(["region", "feature_type"]).size())
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
