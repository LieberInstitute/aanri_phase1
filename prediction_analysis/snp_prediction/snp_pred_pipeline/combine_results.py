#!/usr/bin/env python

import pandas as pd
from glob import iglob


def map_tissues(tissue):
    return {'caudate': "Caudate", "dlpfc": "DLPFC",
            "hippocampus": "Hippocampus", "dentateGyrus": "Dentate Gyrus"}[tissue]


def load_data(fn, tissue):
    df = pd.read_csv("../../_m/%s/%s" % (tissue, fn), sep='\t')
    df["Tissue"] = map_tissues(tissue)
    return df


def partial_R2_metrics():
    cc = load_data("raffe_partial_r2.tsv", "caudate")
    dd = load_data("raffe_partial_r2.tsv", "dlpfc")
    hh = load_data("raffe_partial_r2.tsv", "hippocampus")
    gg = load_data("raffe_partial_r2.tsv", "dentateGyrus")
    pd.concat([cc, dd, hh, gg], axis=0)\
      .to_csv("rf_partial_r2_metrics.tsv", sep='\t', index=False)


def partial_R2_metrics_enet():
    cc = load_data("enet_partial_r2.tsv", "caudate")
    dd = load_data("enet_partial_r2.tsv", "dlpfc")
    hh = load_data("enet_partial_r2.tsv", "hippocampus")
    gg = load_data("enet_partial_r2.tsv", "dentateGyrus")
    pd.concat([cc, dd, hh, gg], axis=0)\
      .to_csv("enet_partial_r2_metrics.tsv", sep='\t', index=False)


def combine_results(fn, tissue):
    df = pd.DataFrame()
    for filename in iglob("../../_m/%s/*/%s" % (tissue, fn)):
        tmp_df = pd.read_csv(filename, sep='\t')
        df = pd.concat([df, tmp_df], axis=0)
    return df


def idv_partial():
    cc = combine_results("individual_partial_r2.tsv", "caudate")
    dd = combine_results("individual_partial_r2.tsv", "dlpfc")
    hh = combine_results("individual_partial_r2.tsv", "hippocampus")
    gg = combine_results("individual_partial_r2.tsv", "dentateGyrus")
    pd.concat([cc, dd, hh, gg], axis=0)\
      .to_csv("individual_partial_r2_metrics.tsv", sep='\t', index=False)


def main():
    partial_R2_metrics()
    partial_R2_metrics_enet()
    idv_partial()    


if __name__ == '__main__':
    main()
