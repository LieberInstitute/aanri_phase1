"""
This script runs tensorQTL in python.
"""
import pandas as pd
import argparse, genotypeio
from functools import lru_cache
from tensorqtl import susie, read_phenotype_bed

@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../plink_format/_m/genotypes"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(fam_id="fid"), variant_df


@lru_cache()
def get_interaction():
    interaction_file = "../../_m/ancestry_interaction.txt"
    return pd.read_csv(interaction_file, sep='\t', index_col=0).loc[:,"Eur"]


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

    # Run SuSiE fine mapping
    susie_df = susie.map(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                         covariates_df=get_covars(args.feature), max_iter=200,
                         maf_threshold=0.05, window=500000)
    susie_df.to_csv("%s.susie_finemapped.txt.gz" % prefix, sep='\t', index=False)


if __name__ == "__main__":
    main()
