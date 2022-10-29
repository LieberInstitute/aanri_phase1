"""
This script summarized DE results using mash model.
"""
import session_info
import pandas as pd
from functools import lru_cache

@lru_cache()
def get_mash_eqtl(feature, tissue, fdr):
    # lfsr < 0.05 for significant eQTL by features
    df = pd.read_csv("../../%s/_m/lfsr_allpairs_ancestry.txt.gz"% feature, sep='\t')\
           .loc[:, ["effect", "gene_id", "variant_id", tissue]]\
           .rename(columns={tissue: "lfsr"})
    return df[(df["lfsr"] < fdr)].copy()


@lru_cache()
def get_mash_es(feature, tissue):
    # get effect size
    return pd.read_csv("../../%s/_m/posterior_mean_allpairs_ancestry.txt.gz"%feature,
                     sep='\t')\
             .loc[:, ["gene_id", "variant_id", tissue]]\
             .rename(columns={tissue: "posterior_mean"})


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
def annotate_eqtl(feature, tissue, fdr):
    return get_mash_eqtl(feature, tissue, fdr)\
        .merge(get_mash_es(feature, tissue), on=["gene_id", "variant_id"])\
        .merge(get_annotation(feature), left_on="gene_id", right_on="names")\
        .drop(["names"], axis=1)


@lru_cache()
def extract_features(tissue, fdr):
    # Extract DE from mash model
    genes = annotate_eqtl("genes", tissue, fdr)
    trans = annotate_eqtl("transcripts", tissue, fdr)
    exons = annotate_eqtl("exons", tissue, fdr)
    juncs = annotate_eqtl("junctions", tissue, fdr)
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
        for variable in ["effect", "gene_id", "gencodeID"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print("\nGene:\t\t%d\nTranscript:\t%d\nExon:\t\t%d\nJunction:\t%d\n" %
                  (gg, tt, ee, jj), file=f)


def get_eqtl_result_by_tissue(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    genes["Feature"] = "Gene"
    trans["Feature"] = "Transcript"
    exons["Feature"] = "Exon"
    juncs["Feature"] = "Junction"
    df = pd.concat([genes, trans, exons, juncs])
    df["Feature"] = df.Feature.astype("category")\
                   .cat.reorder_categories(["Gene", "Transcript", "Exon", "Junction"])
    df["Tissue"] = tissue.replace("_", " ")
    return df


def main():
    bigdata1 = []
    for tissue in ["Caudate", "Dentate_Gyrus", "DLPFC", "Hippocampus"]:
        print_summary(tissue)
        data1 = get_eqtl_result_by_tissue(tissue)
        bigdata1.append(data1)
    df1 = pd.concat(bigdata1)
    cols = ["Tissue", "gene_id", "variant_id", "gencodeID", "Symbol",
            "seqnames", "start", "end", "lfsr", "posterior_mean", "Feature"]
    df1.sort_values(["Tissue", "Feature", "lfsr", "posterior_mean"])\
       .loc[:, cols]\
       .to_csv("BrainSeq_ancestry_4features_4regions.txt.gz",
               sep='\t', index=False)
    ## Reproducibilty information
    session_info.show()


if __name__ == "__main__":
    main()
