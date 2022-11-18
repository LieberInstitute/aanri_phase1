"""
This script runs tensorQTL in python.
"""
import pandas as pd
from pyhere import here
from functools import lru_cache
import argparse, torch, session_info
from tensorqtl import read_phenotype_bed, genotypeio, Residualizer

def map_tissue(tissue):
    return {"caudate": "Caudate", "dlpfc": "DLPFC",
            "dentateGyrus": "Dentate_Gyrus",
            "hippocampus": "Hippocampus"}[tissue]


def map_feature(feature):
    return {"genes": "Gene", "transcripts": "Transcript",
            "exons": "Exon", "junctions": "Junction"}[feature]


def get_phenotype(tissue, feature):
    expr_bed = here(f"eqtl_analysis/{tissue}/{feature}/",
                    "normalize_expression/_m/",
                    f"{feature}.expression.bed.gz")
    return read_phenotype_bed(str(expr_bed))


@lru_cache()
def get_residualized(tissue, feature):
    # Load residualized data
    s_file = here("input/phenotypes/merged/_m/merged_phenotypes.csv")
    samp   = pd.read_csv(s_file, index_col=0)\
               .loc[:, ["BrNum"]]
    res_file = here(f"differential_analysis/{tissue}/_m/",
                    f"{feature}/residualized_expression.tsv")
    return pd.read_csv(res_file, index_col=0, sep='\t')\
             .rename(columns=samp.to_dict()["BrNum"])


@lru_cache()
def load_de(tissue, feature):
    de_file = here("differential_analysis/tissue_comparison/summary_table/_m/",
                   "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    de_df = pd.read_csv(de_file, sep='\t')
    return de_df[(de_df["Tissue"] == map_tissue(tissue).replace("_", " ")) &
                 (de_df["Type"] == map_feature(feature))].copy()


@lru_cache()
def subset_phenotypes(tissue, feature):
    phenotype_df, phenotype_pos_df = get_phenotype(tissue, feature)
    res_df = get_residualized(tissue, feature)
    shared_samples = set(res_df.columns) & \
        set(phenotype_df.columns)
    shared_features = set(res_df.index) & \
        set(phenotype_df.index)
    return res_df.loc[shared_features, shared_samples].sort_index()


@lru_cache()
def get_effectsize(tissue, feature):
    es_file = here(f"eqtl_analysis/tissue_comparison/{feature}/_m/",
                   "posterior_mean_allpairs_ancestry.txt.gz")
    return pd.read_csv(es_file, sep='\t')\
           .loc[:, ["effect", "gene_id", "variant_id", map_tissue(tissue)]]


@lru_cache()
def load_eFeatures(tissue, feature):
    egene_file = here(f"eqtl_analysis/{tissue}/{feature}/",
                      "cis_analysis/_m/LIBD_TOPMed_AA.genes.txt.gz")
    cols = ["phenotype_id","variant_id","slope","slope_se",
            "pval_nominal","qval"]
    return pd.read_csv(egene_file, sep='\t').loc[:, cols]


@lru_cache()
def get_genotypes(tissue, feature):
    plink_prefix_path = here(f"eqtl_analysis/{tissue}/{feature}/",
                             "plink_format/_m/genotypes")
    pr = genotypeio.PlinkReader(str(plink_prefix_path))
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(fam_id="fid"), variant_df


@lru_cache()
def subset_genotypes(tissue, feature):
    genotype_df, variant_df = get_genotypes(tissue, feature)
    return genotype_df.loc[load_eFeatures(tissue, feature).variant_id]


@lru_cache()
def get_topQTL(tissue, feature):
    es_df        = get_effectsize(tissue, feature)
    genotype_df  = subset_genotypes(tissue, feature)
    phenotype_df = subset_phenotypes(tissue, feature)
    eFeatures    = load_eFeatures(tissue, feature)\
        .loc[:, ["phenotype_id", "variant_id", "qval"]]\
        .merge(es_df, left_on=["phenotype_id", "variant_id"],
               right_on=["gene_id", "variant_id"])\
        .drop(["gene_id", "effect"], axis=1)\
        .rename(columns={map_tissue(tissue): "posterior_mean"})\
        .set_index("phenotype_id")
    shared_features = list(set(phenotype_df.index) & set(eFeatures.index))
    top_feature = eFeatures.sort_values("qval")\
                           .loc[shared_features, ["variant_id","qval",
                                                  "posterior_mean"]]
    return pd.merge(top_feature, genotype_df, left_on="variant_id",
                    right_index=True).reset_index().drop_duplicates()


def get_predExpr_topQTL(tissue, feature):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    beta = get_topQTL(tissue, feature).set_index("phenotype_id")\
                                      .loc[:, ["posterior_mean"]]
    top_df = get_topQTL(tissue, feature)\
        .set_index("phenotype_id")\
        .drop(["variant_id","qval", "posterior_mean"], axis=1)
    genotypes_df = torch.tensor(top_df.values,
                                dtype=torch.float32).to(device)
    beta_df = torch.tensor(beta.values,
                           dtype=torch.float32).to(device)
    # Predict expression
    return pd.DataFrame(torch.mul(beta_df, genotypes_df).cpu(),
                        columns=top_df.columns,
                        index=top_df.index).astype("float")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tissue', type=str)
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()

    # Extract predicted expression (Effect Size X Dosage [0|1|2])
    tissue = args.tissue; feature = args.feature
    get_predExpr_topQTL(tissue, feature)\
        .to_csv(f"{feature}.predicted_expression.{tissue}.txt.gz", sep='\t')
    ## Session information
    session_info.show()
    

if __name__ == "__main__":
    main()
