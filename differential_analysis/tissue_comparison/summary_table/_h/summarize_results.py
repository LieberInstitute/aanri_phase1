"""
This script summarized DE results using mash model.
"""
import pandas as pd
from functools import lru_cache


@lru_cache()
def get_mash_degs(feature, tissue):
    # lfsr < 0.05 for significant DE features
    df = pd.read_csv("../../_m/%s/lfsr_feature_4tissues.txt.gz"% feature, sep='\t')\
           .loc[:, ["Effect", tissue]]
    return df[(df[tissue] < 0.05)].rename(columns={tissue: "lfsr"})


@lru_cache()
def get_mash_es(feature, tissue):
    # get effect size
    return pd.read_csv("../../_m/%s/posterior_mean_feature_4tissues.txt.gz"%feature,
                     sep='\t')\
             .loc[:, ["Effect",tissue]].rename(columns={tissue: "posterior_mean"})


@lru_cache()
def get_annotation(feature):
    base_loc = "/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/input/text_files_counts/"
    config = {
        "genes": "%s/_m/caudate/gene_annotation.tsv" % base_loc,
        "transcripts": "%s/_m/caudate/tx_annotation.tsv" % base_loc,
        "exons": "%s/_m/caudate/exon_annotation.tsv" % base_loc,
        "junctions": "%s/_m/caudate/jxn_annotation.tsv" % base_loc,
    }
    return pd.read_csv(config[feature], sep='\t')\
             .loc[:, ["names", "seqnames", "start", "end", "Symbol", "gencodeID"]]


@lru_cache()
def annotate_degs(feature, tissue):
    return get_mash_degs(feature, tissue)\
        .merge(get_mash_es(feature, tissue), on="Effect")\
        .merge(get_annotation(feature), left_on="Effect", right_on="names")\
        .drop(["names"], axis=1)


@lru_cache()
def extract_features(tissue):
    # Extract DE from mash model
    genes = annotate_degs("genes", tissue)
    trans = annotate_degs("transcripts", tissue)
    exons = annotate_degs("exons", tissue)
    juncs = annotate_degs("junctions", tissue)
    return genes, trans, exons, juncs


def print_summary(tissue):
    genes, trans, exons, juncs = extract_features(tissue)
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


def get_DEGs_result_by_tissue(tissue):
    genes, trans, exons, juncs = extract_features(tissue)
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
    bigdata = []
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        print_summary(tissue)
        data = get_DEGs_result_by_tissue(tissue)
        bigdata.append(data)
    df = pd.concat(bigdata)
    cols = ["Tissue", "Effect", "gencodeID", "Symbol", "seqnames", "start", "end",
            "lfsr", "posterior_mean", "Type"]
    df.sort_values(["Tissue", "Type", "lfsr", "posterior_mean"]).loc[:, cols]\
      .to_csv("BrainSeq_ancestry_4features_4regions.txt.gz",
              sep='\t', index=False)


if __name__ == "__main__":
    main()
