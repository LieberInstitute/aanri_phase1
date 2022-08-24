"""
This script runs tensorQTL in python.
"""
import pandas as pd
import argparse, genotypeio
from functools import lru_cache
from tensorqtl import cis, read_phenotype_bed, calculate_qvalues

@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../plink_format/_m/genotypes"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(fam_id="fid"), variant_df


@lru_cache()
def get_covars(feature):
    covar_file = "../../covariates/_m/%s.combined_covariates.txt" % feature
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


@lru_cache()
def get_phenotype(feature):
    expr_bed = "../../normalize_expression/_m/%s.expression.bed.gz" % feature
    return read_phenotype_bed(expr_bed)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()

    # Load data
    phenotype_df, phenotype_pos_df = get_phenotype(args.feature)
    genotype_df, variant_df = get_genotypes()
    prefix = "LIBD_TOPMed_AA"

    # Nominal
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                    prefix, covariates_df=get_covars(args.feature),
                    maf_threshold=0.01, window=500000, output_dir=".")
    # Permutation
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df,phenotype_pos_df,
                         covariates_df=get_covars(args.feature),
                         maf_threshold=0.01, window=500000, seed=13131313)
    calculate_qvalues(cis_df, fdr=0.05)
    cis_df.to_csv("%s.genes.txt.gz" % prefix, sep='\t')
    # Conditional
    indep_df = cis.map_independent(genotype_df, variant_df, cis_df,phenotype_df,
                                   phenotype_pos_df, covariates_df=get_covars(args.feature),
                                   maf_threshold=0.01, window=500000, seed=13131313)
    indep_df.to_csv("%s.conditional.txt.gz" % prefix, sep='\t')


if __name__ == "__main__":
    main()
