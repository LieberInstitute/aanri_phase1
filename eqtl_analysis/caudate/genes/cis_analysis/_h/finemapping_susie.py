"""
This script runs tensorQTL in python.
"""
import pandas as pd
from functools import lru_cache
import argparse, torch, genotypeio
from tensorqtl import susie, read_phenotype_bed
from tensorqtl.post import get_significant_pairs

print(f"PyTorch {torch.__version__}")
print(f"Pandas {pd.__version__}")

@lru_cache()
def get_covars(feature):
    covar_file = "../../covariates/_m/%s.combined_covariates.txt" % feature
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


@lru_cache()
def get_phenotype(feature):
    expr_bed = "../../normalize_expression/_m/%s.expression.bed.gz" % feature
    return read_phenotype_bed(expr_bed)


@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../plink_format/_m/genotypes"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(fam_id="fid"), variant_df


@lru_cache()
def get_permutation_results(prefix):
    return pd.read_csv("%s.genes.txt.gz" % (prefix), sep="\t", index_col=0)


@lru_cache()
def get_eGenes(prefix):
    df = get_permutation_results(prefix)\
        .loc[:, ["qval", "pval_nominal_threshold"]]
    return df[(df["qval"] <= 0.05)].copy()


@lru_cache()
def get_eqtl(prefix):
    return get_significant_pairs(get_permutation_results(prefix), prefix)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str, default="genes")
    parser.add_argument('--ancestry', type=str, default="AA")
    args=parser.parse_args()

    # Load data
    phenotype_df, phenotype_pos_df = get_phenotype(args.feature)
    genotype_df, variant_df = get_genotypes()
    prefix = "LIBD_TOPMed_%s" % args.ancestry
    covariates_df = get_covars(args.feature)

    # Load eQTL results
    egene_df = get_eGenes(prefix).reset_index()
    eqtl_df = get_eqtl(prefix)

    # Filter for cis-eQTL (eFeatures)
    phenotype_df = phenotype_df.loc[egene_df.phenotype_id,:]
    phenotype_pos_df = phenotype_pos_df.loc[egene_df.phenotype_id,:]
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
