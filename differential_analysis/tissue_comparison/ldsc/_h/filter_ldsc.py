## This script reviews LDSC results.

import pandas as pd
from functools import lru_cache
from statsmodels.stats.multitest import fdrcorrection

@lru_cache()
def load_ldsc_results():
    return pd.read_excel("ldsc_results_final.xlsx")


@lru_cache()
def preprocess_data():
    df = load_ldsc_results()
    group_df = df.group\
                 .str\
                 .split("_", expand=True)\
                 .rename(columns={0:"Tissue",1:"Feature",2:"Direction"})
    return pd.concat([group_df, df], axis=1)\
             .drop(["group"], axis=1)


def cal_fdr(tissue, feature):
    df0 = df = preprocess_data()
    df = df0[(df0["Feature"] == feature) &
             (df0["Tissue"] == tissue)].copy()
    _, fdr = fdrcorrection(df.Enrichment_p)
    df["Enrichment_FDR"] = fdr
    return df[(df["Prop._SNPs"] > 0.01)].copy()


def filter_data():
    new_df = pd.DataFrame()
    for tissue in ["Caudate", "DentateGyrus", "DLPFC", "Hippocampus"]:
        for feature in ["Gene", "Transcript"]:
            new_df = pd.concat([new_df,
                                cal_fdr(tissue, feature)], axis=0)
    return new_df


def split_feature(df):
    gg = df[(df["Feature"] == "Gene")].copy()
    tt = df[(df["Feature"] == "Transcript")].copy()
    return gg, tt


def split_regions(df):
    cc = df[(df["Tissue"] == "Caudate")].copy()
    gg = df[(df["Tissue"] == "DentateGyrus")].copy()
    dd = df[(df["Tissue"] == "DLPFC")].copy()
    hh = df[(df["Tissue"] == "Hippocampus")].copy()
    return cc, gg, dd, hh


def shared_traits(df):
    cc, gg, dd, hh = split_regions(df)
    cc_t = cc.loc[(cc["Enrichment_FDR"] < 0.05), "trait"]
    gg_t = gg.loc[(gg["Enrichment_FDR"] < 0.05), "trait"]
    dd_t = dd.loc[(dd["Enrichment_FDR"] < 0.05), "trait"]
    hh_t = hh.loc[(hh["Enrichment_FDR"] < 0.05), "trait"]
    return list(set(cc_t) | set(gg_t) | set(dd_t) | set(hh_t))


def generate_data():
    gene, tx = split_feature(filter_data())
    sig_traits = shared_traits(tx)
    gene = gene[(gene["trait"].isin(sig_traits))].copy()
    tx = tx[(tx["trait"].isin(sig_traits))].copy()
    return pd.concat([gene,tx], axis=0)


def main():
    generate_data()\
        .to_csv("filtered_ldsc_results.tsv", sep='\t', index=False)


if __name__ == "__main__":
    main()
