"""
This script generates venn diagrams across brain region for each
direction and dataset for both DEG and TWAS.
"""
import functools
import numpy as np
import pandas as pd
from venn import venn
import matplotlib.pyplot as plt

@functools.lru_cache()
def get_data():
    return pd.read_csv("../../_m/clinical_overlap_ancestryDEGs.txt.gz",
                       sep='\t')


def subset_data(dataset, direction, method):
    df = get_data()[(get_data()["Direction"]==direction) &
                    (get_data()["Method"]==method) &
                    (get_data()["Dataset"]==dataset)]
    cc = set(df[(df["Tissue"]=="Caudate")].Genes)
    gg = set(df[(df["Tissue"]=="Dentate Gyrus")].Genes)
    dd = set(df[(df["Tissue"]=="DLPFC")].Genes)
    hh = set(df[(df["Tissue"]=="Hippocampus")].Genes)
    return cc, gg, dd, hh


def plot_venn(dataset, direction, method):
    cc, gg, dd, hh = subset_data(dataset, direction, method)
    tissues = {
        "Caudate":cc, "Dentate Gyrus":gg,
        "DLPFC": dd, "Hippocampus": hh
    }
    fn = "%s_%s_%s.pdf" % (dataset, direction, method)
    venn(tissues, fmt="{size}", fontsize=15, legend_loc="lower right")
    plt.title(dataset.replace("_", " "), fontweight="bold", size=20)
    plt.savefig(fn.replace(" ", "_"))


def main():
    for method in get_data().Method.unique():
        for direction in get_data().Direction.unique():
            for dataset in get_data().Dataset.unique():
                try:
                    plot_venn(dataset, direction, method)
                except:
                    None


if __name__ == "__main__":
    main()
