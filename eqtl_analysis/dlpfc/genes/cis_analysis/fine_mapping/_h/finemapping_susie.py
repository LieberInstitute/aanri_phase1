"""
This script runs tensorQTL in python.
"""
import pandas as pd
import argparse, torch
from pyhere import here
from functools import lru_cache
from tensorqtl.post import get_significant_pairs
from tensorqtl import susie, read_phenotype_bed, genotypeio

print(f"PyTorch {torch.__version__}")
print(f"Pandas {pd.__version__}")

@lru_cache()
def get_covars(feature):
    covar_file = f"../../../covariates/_m/{feature}.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


@lru_cache()
def get_phenotype(feature):
    expr_bed = f"../../../normalize_expression/_m/{feature}.expression.bed.gz"
    return read_phenotype_bed(expr_bed)


@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../../plink_format/_m/genotypes"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(fam_id="fid"), variant_df


@lru_cache()
def get_permutation_results():
    efile = here("eqtl_analysis/tissue_comparison/feature_summary/_m/",
                 "BrainSeq_ancestry_4features_4regions.txt.gz")
    return pd.read_csv(efile, sep="\t")


@lru_cache()
def get_eqtl():
    df = get_permutation_results()\
        .loc[:, ["Tissue", "gene_id", "variant_id", "Feature"]]
    return df[(df["Tissue"]  == "DLPFC") &
              (df["Feature"] == "Gene")].copy()


@lru_cache()
def get_eGenes():
    return get_eqtl().drop_duplicates(subset="gene_id")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str, default="genes")
    parser.add_argument('--ancestry', type=str, default="AA")
    args=parser.parse_args()

    # Load data
    phenotype_df, phenotype_pos_df = get_phenotype(args.feature)
    genotype_df, variant_df = get_genotypes()
    prefix = f"LIBD_TOPMed_{args.ancestry}"
    covariates_df = get_covars(args.feature)

    # Load eQTL results
    egene_df = get_eGenes()
    eqtl_df = get_eqtl()

    # Filter for cis-eQTL (eFeatures)
    phenotype_df = phenotype_df.loc[egene_df.gene_id,:]
    phenotype_pos_df = phenotype_pos_df.loc[egene_df.gene_id,:]
    genotype_df = genotype_df.loc[eqtl_df.variant_id.unique(), :]
    variant_df = variant_df.loc[eqtl_df.variant_id.unique(), :]

    # Run SuSiE fine mapping by chromosome
    chroms = sorted(phenotype_pos_df.chr.unique())
    susie_df = []
    for k,chrom in enumerate(chroms, 1):
        print(f'  * processing chr. {k}/{len(chroms)}', flush=True)
        dfx = susie.map(genotype_df, variant_df,
                        phenotype_df.loc[phenotype_pos_df['chr']==chrom],
                        phenotype_pos_df.loc[phenotype_pos_df['chr']==chrom],
                        covariates_df=covariates_df, max_iter=200,
                        maf_threshold=0.01, window=500000)
        susie_df.append(dfx)
    pd.concat(susie_df, axis=0)\
      .to_csv("%s.susie.txt.gz" % prefix, sep='\t', index=False)


if __name__ == "__main__":
    main()
