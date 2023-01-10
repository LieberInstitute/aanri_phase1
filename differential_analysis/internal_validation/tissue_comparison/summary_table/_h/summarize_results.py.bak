"""
This script summarized DE results using mash model.
"""
import pandas as pd
from functools import lru_cache


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
def annotate_degs(feature, tissue, fdr):
    return get_mash_degs(feature, tissue, fdr)\
        .merge(get_mash_es(feature, tissue), on="Effect")\
        .merge(get_annotation(feature), left_on="Effect", right_on="names")\
        .drop(["names"], axis=1)


@lru_cache()
def extract_features(tissue, fdr):
    # Extract DE from mash model
    genes = annotate_degs("genes", tissue, fdr)
    trans = annotate_degs("transcripts", tissue, fdr)
    exons = annotate_degs("exons", tissue, fdr)
    juncs = annotate_degs("junctions", tissue, fdr)
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
    cols = ["Tissue", "Effect", "gencodeID", "Symbol", "seqnames", "start", "end",
            "lfsr", "posterior_mean", "Type"]
    df1.sort_values(["Tissue", "Type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions.txt.gz",sep='\t', index=False)
    df2.sort_values(["Tissue", "Type", "lfsr", "posterior_mean"]).loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz",
               sep='\t', index=False)


if __name__ == "__main__":
    main()
