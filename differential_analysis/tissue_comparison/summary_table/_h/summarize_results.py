"""
This script summarized DE results using mash model.
"""
import pandas as pd
import session_info
from pyhere import here
from gtfparse import read_gtf
from functools import lru_cache

@lru_cache()
def get_gtf():
    gtf_file = "/dcl02/lieber/apaquola/genome/human/gencode25/gtf.CHR/"+\
        "_m/gencode.v25.annotation.gtf"
    return read_gtf(gtf_file)


def gene_annotation():
    gtf = get_gtf()[get_gtf()["feature"] == "transcript"]
    return gtf[["transcript_id", "gene_id", "gene_name", "gene_type",
                "seqname", "start", "end", "strand"]]

@lru_cache()
def get_mash_degs(feature, tissue, fdr):
    # lfsr < 0.05 for significant DE features
    df = pd.read_csv("../../_m/%s/lfsr_feature_4tissues.txt.gz"% feature, sep='\t')\
           .loc[:, ["Effect", tissue]]
    return df[(df[tissue] < fdr)].rename(columns={tissue: "lfsr"})


@lru_cache()
def get_mash_es(feature, tissue):
    # get effect size
    return pd.read_csv("../../_m/%s/posterior_mean_feature_4tissues.txt.gz"%feature,
                     sep='\t')\
             .loc[:, ["Effect",tissue]].rename(columns={tissue: "posterior_mean"})


@lru_cache()
def get_annotation(feature):
    config = {
        "genes": here("input/text_files_counts/_m",
                      "caudate/gene_annotation.tsv"),
        "transcripts": here("input/text_files_counts/_m",
                            "caudate/tx_annotation.tsv"),
        "exons": here("input/text_files_counts/_m",
                      "caudate/exon_annotation.tsv"),
        "junctions": here("input/text_files_counts/_m",
                          "caudate/jxn_annotation.tsv"),
    }
    return pd.read_csv(config[feature], sep='\t')\
             .loc[:, ["names", "seqnames", "start", "end", "Symbol", "gencodeID"]]


@lru_cache()
def annotate_degs(feature, tissue, fdr):
    return get_mash_degs(feature, tissue, fdr)\
        .merge(get_mash_es(feature, tissue), on="Effect")\
        .merge(get_annotation(feature), left_on="Effect", right_on="names")\
        .drop(["names"], axis=1)


@lru_cache()
def extract_features(tissue, fdr):
    # Gene annotation
    gtf_annot = gene_annotation()
    annot_gene = gtf_annot.loc[:, ["gene_id", "gene_type", "strand"]].drop_duplicates()
    annot_tx   = gtf_annot.loc[:, ["transcript_id", "gene_type", "strand"]]
    # Extract DE from mash model
    genes = annotate_degs("genes", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeID", right_on="gene_id", how="left")\
        .drop("gene_id", axis=1)
    trans = annotate_degs("transcripts", tissue, fdr)\
        .merge(annot_tx, left_on="Effect", right_on="transcript_id", how="left")\
        .drop("transcript_id", axis=1)
    exons = annotate_degs("exons", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeID", right_on="gene_id", how="left")\
        .drop("gene_id", axis=1)
    juncs = annotate_degs("junctions", tissue, fdr)\
        .merge(annot_gene, left_on="gencodeID", right_on="gene_id", how="left")\
        .drop("gene_id", axis=1)
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
        for variable in ["Effect", "gencodeID"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print("\nGene:\t\t%d\nTranscript:\t%d\nExon:\t\t%d\nJunction:\t%d\n" %
                  (gg, tt, ee, jj), file=f)


def get_DEGs_result_by_tissue(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    genes["Type"] = "Gene"
    trans["Type"] = "Transcript"
    exons["Type"] = "Exon"
    juncs["Type"] = "Junction"
    df = pd.concat([genes, trans, exons, juncs])
    df["Type"] = df.Type.astype("category")\
                   .cat.reorder_categories(["Gene", "Transcript", "Exon", "Junction"])
    df["Tissue"] = tissue.replace(".", " ")
    return df


def main():
    bigdata1 = []; bigdata2 = [];
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        print_summary(tissue)
        data1 = get_DEGs_result_by_tissue(tissue)
        data2 = get_DEGs_result_by_tissue(tissue, 1)
        bigdata1.append(data1); bigdata2.append(data2)
    df1 = pd.concat(bigdata1); df2 = pd.concat(bigdata2)
    # Summary
    with open("effect_sizes.log", mode="w") as f:
        print("Effect size:", file=f)
        print(df1.loc[:, ["Tissue", "Type", "posterior_mean"]]\
              .groupby(["Tissue", "Type"]).describe().to_string(), file=f)
    print("\nSummary:")
    gene = df1[(df1["Type"] == "Gene")].drop_duplicates(subset="gencodeID")
    print(gene.shape)
    print(gene.groupby("gene_type").size())
    # Output
    cols = ["Tissue", "Effect", "gencodeID", "Symbol", "seqnames", "start", "end",
            "strand", "gene_type", "lfsr", "posterior_mean", "Type"]
    df1.sort_values(["Tissue", "Type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions.txt.gz",sep='\t', index=False)
    df2.sort_values(["Tissue", "Type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz",
               sep='\t', index=False)
    # Session infomation
    session_info.show()


if __name__ == "__main__":
    main()
