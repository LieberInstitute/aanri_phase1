"""
This script summarized DE results using mash model.
"""
import pandas as pd
import polars as pl
import session_info
from pyhere import here
from gtfparse import read_gtf
from functools import lru_cache

@lru_cache()
def get_annotation(feature):
    config = {
        "genes": here("input/text_files_counts/_m/caudate/gene_annotation.tsv"),
        "transcripts": here("input/text_files_counts/_m/caudate/tx_annotation.tsv"),
        "exons": here("input/text_files_counts/_m/caudate/exon_annotation.tsv"),
        "junctions": here("input/text_files_counts/_m/caudate/jxn_annotation.tsv"),
    }
    return pd.read_csv(config[feature], sep='\t')\
             .loc[:, ["names", "seqnames", "start", "end", "Symbol", "gencodeID"]]\
             .rename(columns=str.lower)


@lru_cache()
def get_gtf():
    gtf_file = "/dcl02/lieber/apaquola/genome/human/gencode25/gtf.CHR/"+\
        "_m/gencode.v25.annotation.gtf"
    return read_gtf(gtf_file)


def gene_annotation():
    gtf = get_gtf().filter(pl.col("feature") == "transcript")
    return gtf[["transcript_id", "gene_id", "gene_name", "gene_type",
                "seqname", "start", "end", "strand"]]

@lru_cache()
def get_mash_degs(feature, tissue, fdr):
    if tissue == "Dentate Gyrus":
        new_tissue = "dentateGyrus"
    else:
        new_tissue = tissue.lower()
    df = pd.read_csv(f"../../_m/{feature}/mash_lfsr_local.txt", sep='\t')\
           .loc[:, ["feature_id", new_tissue]]
    return df[(df[new_tissue] < fdr)].rename(columns={new_tissue: "lfsr"})


@lru_cache()
def get_mash_es(feature, tissue):
    if tissue == "Dentate Gyrus":
        new_tissue = "dentateGyrus"
    else:
        new_tissue = tissue.lower()
    df = pd.read_csv(f"../../_m/{feature}/mash_effectsize_local.txt",
                     sep='\t')\
           .loc[:, ["feature_id", new_tissue]]\
           .rename(columns={new_tissue: "posterior_mean"})
    ## Flip sign to match global ancestry
    df["posterior"] = df["posterior_mean"] * -1
    return df


@lru_cache()
def annotate_degs(feature, tissue, fdr):
    return get_mash_degs(feature, tissue, fdr)\
        .merge(get_mash_es(feature, tissue), on="feature_id")\
        .merge(get_annotation(feature), left_on="feature_id", right_on="names")\
        .drop(["names"], axis=1)


@lru_cache()
def extract_features(tissue, fdr):
    # Gene annotation
    gtf_annot = gene_annotation()
    annot_gene = gtf_annot[:, ["gene_id", "gene_type"]]\
                          .unique().to_arrow().to_pandas()
    annot_tx   = gtf_annot[:, ["transcript_id", "gene_id", "gene_type"]]\
        .to_arrow().to_pandas()
    # Extract DE from mash model
    genes = annotate_degs("genes", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeid", right_on="gene_id", how="left")\
        .drop("gencodeid", axis=1)
    trans = annotate_degs("transcripts", tissue, fdr)\
        .merge(annot_tx, left_on="feature_id", right_on="transcript_id", how="left")\
        .drop(["gencodeid", "transcript_id"], axis=1)
    exons = annotate_degs("exons", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeid", right_on="gene_id", how="left")\
        .drop("gencodeid", axis=1)
    juncs = annotate_degs("junctions", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeid", right_on="gene_id", how="left")\
        .drop("gencodeid", axis=1)
    return genes, trans, exons, juncs


def print_summary(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    if tissue == "Caudate":
        w_mode = "w"
    else:
        w_mode = "a"
    statement = "Significant DE (lfsr < 0.05) in %s" % tissue
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for variable in ["feature_id", "gene_id"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print("\nGene:\t\t%d\nTranscript:\t%d\nExon:\t\t%d\nJunction:\t%d\n" %
                  (gg, tt, ee, jj), file=f)


def get_DEGs_result_by_tissue(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    genes["feature_type"] = "Gene"
    trans["feature_type"] = "Transcript"
    exons["feature_type"] = "Exon"
    juncs["feature_type"] = "Junction"
    df = pd.concat([genes, trans, exons, juncs])
    df["feature_type"] = df.feature_type\
                           .astype("category")\
                           .cat.reorder_categories(["Gene", "Transcript",
                                                    "Exon", "Junction"])
    df["region"] = tissue
    return df


@lru_cache()
def extract_partial(tissue, fdr=0.05):
    gtf_annot = gene_annotation()
    annot_gene = gtf_annot[:, ["gene_id", "gene_type"]]\
                          .unique().to_arrow().to_pandas()
    annot_tx   = gtf_annot[:, ["transcript_id", "gene_id", "gene_type"]]\
        .to_arrow().to_pandas()
    # Extract DE from mash model
    genes = annotate_degs("genes", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeid", right_on="gene_id", how="left")\
        .drop("gencodeid", axis=1)
    trans = annotate_degs("transcripts", tissue, fdr)\
        .merge(annot_tx, left_on="feature_id", right_on="transcript_id", how="left")\
        .drop(["gencodeid", "transcript_id"], axis=1)
    exons = annotate_degs("exons", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeid", right_on="gene_id", how="left")\
        .drop("gencodeid", axis=1)
    return genes, trans, exons


def print_partial(tissue, fdr=0.05):
    genes, trans, exons = extract_partial(tissue, fdr)
    if tissue == "Caudate":
        w_mode = "w"
    else:
        w_mode = "a"
    statement = "Significant DE (lfsr < 0.05) in %s" % tissue
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for variable in ["feature_id", "gene_id"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            print(f"\nGene:\t\t{gg}\nTranscript:\t{tt}"+\
                  f"\nExon:\t\t{ee}\n", file=f)


def get_DEGs_result_by_tissue_partial(tissue, fdr=0.05):
    genes, trans, exons = extract_partial(tissue, fdr)
    genes["feature_type"] = "Gene"
    trans["feature_type"] = "Transcript"
    exons["feature_type"] = "Exon"
    df = pd.concat([genes, trans, exons])
    df["feature_type"] = df.feature_type\
                           .astype("category")\
                           .cat.reorder_categories(["Gene", "Transcript",
                                                    "Exon"])
    df["region"] = tissue
    return df


def combine_all_feature():
    bigdata1 = []; bigdata2 = [];
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        print_summary(tissue)
        data1 = get_DEGs_result_by_tissue(tissue)
        data2 = get_DEGs_result_by_tissue(tissue, 1)
        bigdata1.append(data1); bigdata2.append(data2)
    df1 = pd.concat(bigdata1); df2 = pd.concat(bigdata2)
    # Summary
    with open("effect_sizes.log", mode="w") as f:
        print("Effect size:", file=f)
        print(df1.loc[:, ["region", "feature_type", "posterior_mean"]]\
              .groupby(["region", "feature_type"]).describe().to_string(), file=f)
    print("\nSummary:")
    gene = df1[(df1["feature_type"] == "Gene")].drop_duplicates(subset="gene_id")
    print(gene.shape)
    print(gene.groupby("gene_type").size())
    # Output
    cols = ["region", "feature_id", "gene_id", "symbol", "seqnames", "start", "end",
            "lfsr", "posterior_mean", "feature_type"]
    df1.sort_values(["region", "feature_type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions.txt.gz",sep='\t', index=False)
    df2.sort_values(["region", "feature_type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz",
               sep='\t', index=False)


def main():
    # Genes and Transcripts
    bigdata1 = []; bigdata2 = [];
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        print_partial(tissue)
        data1 = get_DEGs_result_by_tissue_partial(tissue)
        data2 = get_DEGs_result_by_tissue_partial(tissue, 1)
        bigdata1.append(data1); bigdata2.append(data2)
    df1 = pd.concat(bigdata1); df2 = pd.concat(bigdata2)
    with open("effect_sizes.log", mode="w") as f:
        print("Effect size:", file=f)
        print(df1.loc[:, ["region", "feature_type", "posterior_mean"]]\
              .groupby(["region", "feature_type"]).describe().to_string(), file=f)
    print("\nSummary:")
    gene = df1[(df1["feature_type"] == "Gene")]\
        .drop_duplicates(subset="gene_id")
    print(gene.shape)
    print(gene.groupby("gene_type").size())
    # Output
    cols = ["region", "feature_id", "gene_id", "symbol", "seqnames", "start", "end",
            "lfsr", "posterior_mean", "feature_type"]
    df1.sort_values(["region", "feature_type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions.txt.gz",sep='\t', index=False)
    df2.sort_values(["region", "feature_type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz",
               sep='\t', index=False)
    # Session infomation
    session_info.show()


if __name__ == "__main__":
    main()
