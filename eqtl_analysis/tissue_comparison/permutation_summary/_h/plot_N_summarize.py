## This script summarizes and plots the significant eQTL features

import pandas as pd
import seaborn as sns
from functools import lru_cache
import matplotlib.pyplot as plt

@lru_cache()
def load_eqtl(tissue, feature):
    filename = "../../../%s/%s/cis_analysis/_m/" % (tissue, feature) +\
        "LIBD_TOPMed_AA.signif_variants.txt.gz"
    df = pd.read_csv(filename, sep='\t', nrows=100,
                     usecols=["phenotype_id","variant_id","pval_nominal",
                              "slope","slope_se"], compression="gzip")
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict(),
                       usecols=["phenotype_id","variant_id","pval_nominal",
                                "slope","slope_se"], compression="gzip")


@lru_cache()
def annotate_eqtl(tissue, feature):
    f_dict = {"genes":"gene", "transcripts":"tx",
              "exons":"exon", "junctions":"jxn"}
    fn = "../../../../input/text_files_counts/_m/caudate/"+\
        "%s_annotation.tsv" % f_dict[feature]
    df = pd.merge(pd.read_csv(fn, sep='\t', index_col=0),
                  load_eqtl(tissue, feature),
                  left_index=True, right_on="phenotype_id")
    cols = ["phenotype_id", "variant_id", "gencodeID", "Symbol", "seqnames",
            "pval_nominal", "slope", "slope_se"]
    return df.loc[:, cols]


@lru_cache()
def combine_features(tissue):
    t_dict = {"caudate": "Caudate", "dentateGyrus": "Dentate_Gyrus",
              "dlpfc": "DLPFC", "hippocampus": "Hippocampus"}
    # Load annotated eQTL
    gene = annotate_eqtl(tissue, "genes"); gene["Feature"] = "Gene"
    tx = annotate_eqtl(tissue, "transcripts"); tx["Feature"] = "Transcript"
    exon = annotate_eqtl(tissue, "exons"); exon["Feature"] = "Exon"
    jxn = annotate_eqtl(tissue, "junctions"); jxn["Feature"] = "Junction"
    # Combine features
    df = pd.concat([gene, tx, exon, jxn], axis=0)
    df["Feature"] = df.Feature.astype("category").cat\
                      .reorder_categories(["Gene", "Transcript",
                                           "Exon", "Junction"])
    df["Region"] = t_dict[tissue]
    return df


@lru_cache()
def merge_regions():
    dt = pd.DataFrame()
    for tissue in ["caudate", "dentateGyrus", "dlpfc", "hippocampus"]:
        dt = pd.concat([dt, combine_features(tissue)], axis=0)
    dt["Region"] = dt.Region.astype("category")
    return dt


def summarize_data():
    dt = merge_regions()
    print("There are a total of %d unique genes across brain regions" % len(dt.gencodeID.unique()))
    print(dt.groupby(["Region", "Feature"]).size())
    print(dt.groupby(["Region", "gencodeID"])\
          .first().reset_index().groupby(["Region", "Feature"]).size())
    print(dt.groupby(["Region", "Symbol"])\
          .first().reset_index().groupby(["Region", "Feature"]).size())
    dt.to_csv("Brainseq_AAonly_signifpairs_4Regions_4Features.txt.gz",
              sep='\t', index=False)


def plotting_features():
    dt = merge_regions().groupby(["Region", "gencodeID"])\
                        .first().reset_index()
    f, ax = plt.subplots(figsize=(6,7))
    bar_plot = sns.countplot(x="Region", data=dt, hue="Feature")
    fig = bar_plot.get_figure()
    fig.savefig("eGenes_across_regions.pdf", bbox_inches="tight")


def main():
    summarize_data()
    plotting_features()


if __name__ == "__main__":
    main()
