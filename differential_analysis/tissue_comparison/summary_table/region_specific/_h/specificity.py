"""
This script extracts the tissue-specific features.
"""
import numpy as np
import pandas as pd
from functools import lru_cache

@lru_cache()
def load_sig():
    return pd.read_csv("../../_m/BrainSeq_ancestry_4features_4regions.txt.gz",
                       sep='\t')


@lru_cache()
def get_tissues(feature):
    cc = load_sig()[(load_sig()["Type"] == feature) &
                    (load_sig()["Tissue"] == "Caudate")].copy()
    dg = load_sig()[(load_sig()["Type"] == feature) &
                    (load_sig()["Tissue"] == "Dentate Gyrus")].copy()
    dd = load_sig()[(load_sig()["Type"] == feature) &
                    (load_sig()["Tissue"] == "DLPFC")].copy()
    hh = load_sig()[(load_sig()["Type"] == feature) &
                    (load_sig()["Tissue"] == "Hippocampus")].copy()
    return cc, dg, dd, hh


@lru_cache()
def region_specific(feature):
    cc, dg, dd, hh = get_tissues(feature)
    caud8 = set(cc.Effect) - set(dg.Effect) - set(dd.Effect) - set(hh.Effect)
    gyrus = set(dg.Effect) - set(cc.Effect) - set(dd.Effect) - set(hh.Effect)
    dlpfc = set(dd.Effect) - set(cc.Effect) - set(dg.Effect) - set(hh.Effect)
    hippo = set(hh.Effect) - set(cc.Effect) - set(dg.Effect) - set(dd.Effect)
    return pd.concat([cc[(cc["Effect"].isin(caud8))],
                      dg[(dg["Effect"].isin(gyrus))],
                      dd[(dd["Effect"].isin(dlpfc))],
                      hh[(hh["Effect"].isin(hippo))]], axis=0)


@lru_cache()
def get_discordant(feature):
    cc, dg, dd, hh = get_tissues(feature)
    shared = set(cc.Effect) & set(dg.Effect) & set(dd.Effect) & set(hh.Effect)
    shared_df = load_sig()[(load_sig()["Effect"].isin(shared))]\
        .loc[:, ["Tissue", "Effect", "posterior_mean"]]\
        .pivot_table(values="posterior_mean",columns=["Tissue"],index=["Effect"])\
        .apply(np.sign)
    discordant = shared_df.eq(shared_df.loc[:, "Caudate"], axis=0).all(axis=1)
    return cc.merge(pd.DataFrame({"Concordant":discordant}).reset_index(),
                    on="Effect")\
             .drop(["Tissue", "lfsr", "posterior_mean"], axis=1)


def main():
    df1 = pd.DataFrame(); df2= pd.DataFrame()
    for feature in ["Gene", "Transcript", "Exon", "Junction"]:
        df1 = pd.concat([df1, region_specific(feature)], axis=0)
        df2 = pd.concat([df2, get_discordant(feature)], axis=0)
    df1.to_csv("BrainSeq_ancestry_AA_region_specific.tsv", sep='\t',index=False)
    df2.to_csv("BrainSeq_ancestry_AA_shared_disconcordant.tsv", sep='\t',
               index=False)


if __name__ == "__main__":
    main()
