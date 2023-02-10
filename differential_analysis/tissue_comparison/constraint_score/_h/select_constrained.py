#!/bin/bash

import pandas as pd

def get_merged_data():
    fn = "brainseq_degs_constrain_score.tsv"
    return pd.read_csv(fn, sep="\t")


def subset_data():
    df = get_merged_data()
    return df[(df["lfsr"] < 0.05)]\
        .sort_values(by=["Tissue", "oe_lof_upper", "lfsr"],
                     ascending=True)


def main():
    df = subset_data()
    cols = ["Tissue", "gencodeID", "gene", "oe_lof_upper_bin",
            "oe_lof_upper", "beta", "lfsr"]
    df.loc[:, cols]\
      .to_csv("BrainSeq_DEG_constrain_score.tsv", sep='\t',
              index=False)


if __name__ == "__main__":
    main()
