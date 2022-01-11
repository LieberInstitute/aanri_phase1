#!/usr/bin/env python
"""
This was developed on Python 3.7+.
"""

import functools
import subprocess
import pandas as pd
import os, argparse, re
from pandas_plink import read_plink

"""
Orginal script generated editted from Apuã. Original
modification by KJ. This script is heavily modified.
"""
__author__ = 'KJ Benjamin'


@functools.lru_cache()
def cached_read_plink(file_prefix, verbose=True):
    return read_plink(file_prefix, verbose)


@functools.lru_cache()
def get_phenotypes(pheno_file):
    """
    Get phenotypes for matching SNPs and expression.
    """
    return pd.read_csv(pheno_file, index_col=0)


@functools.lru_cache()
def get_expression(expr_file):
    """
    Get's the expression for all features. Recommend residualized expression.
    """
    return pd.read_csv(expr_file, sep='\t')


def sample_mapping():
    """
    This fixes the fid numbers for samples.
    """
    fn = "../../_m/Genotype_ID_2222.csv"
    df = pd.read_csv(fn)
    df = df[(df["ID"].str.startswith("R"))].copy()
    df = pd.concat([df, df.ID.str.split("_", expand=True)], axis=1)
    df["IID"] = df[1] + "_" + df[0]
    return df.loc[:, ["BrNum", "IID"]]


def fix_fam(fam):
    """
    Given a fam files from plink, replace the incorrect family ID with
    the correct BrNum.
    """
    ff1 = sample_mapping().merge(fam, left_on="IID", right_on="iid")\
                 .drop(["IID", "fid"], axis=1)\
                 .rename(columns={"BrNum": "fid"})
    ff2 = fam[(fam['fid'] != '0')].copy()
    return pd.concat([ff1, ff2], axis=0)


def select_snps_from_plink_to_df(plink_file_prefix):
    """
    Given a <bim, fam, bed> triple of files, fixes the family IDs (brain number)
    and returns a pandas DataFrame with the SNP genotype values (0,1,2), where
    the rows are family IDs (Br number) and columns are SNP IDs.
    """
    (bim, fam, bed) = cached_read_plink(plink_file_prefix)
    r = pd.DataFrame(bed.compute().transpose())
    new_fam = fix_fam(fam)
    r.index = new_fam['fid']
    r.columns = bim.snp
    return r


def one_hot_encode_snp_df(snp_df):
    """
    Given a snp_df returned by select_snps_from_plink_to_df function, returns
    a one-hot-encoded version of it
    """
    r = pd.get_dummies(snp_df, columns=snp_df.columns, dummy_na=True)
    r.columns = r.columns.str.replace('\.\d+', '', regex=True)
    return r


def add_expression_n_save(expr_file, plink_file_prefix, dirname, pheno_file):
    match = re.search("ENSG\d+.\d+", plink_file_prefix).group(0)
    res_df = pd.merge(get_phenotypes(pheno_file),
                      get_expression(expr_file).loc[match, :],
                      left_index=True, right_index=True)\
               .loc[:, ["BrNum", match]]
    snp_df = select_snps_from_plink_to_df(plink_file_prefix)
    # snp_df.merge(res_df, left_index=True, right_on="BrNum")\
    #       .drop_duplicates(match)\
    #       .drop(["BrNum"], axis=1)\
    #       .to_csv("%s/snps.csv" % dirname)
    oh_snp_df = one_hot_encode_snp_df(snp_df)
    oh_snp_df.merge(res_df, left_index=True, right_on="BrNum")\
             .drop_duplicates(match)\
             .drop(["BrNum"], axis=1)\
             .to_csv("%s/snps_onehot.csv" % dirname)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str)
    parser.add_argument('--expr_file', type=str)
    parser.add_argument('--plink_file_prefix', type=str)
    parser.add_argument('--pheno_file', type=str)
    args=parser.parse_args()
    ## Main section
    os.environ['NUMEXPR_MAX_THREADS'] = '1'
    add_expression_n_save(args.expr_file,
                          args.plink_file_prefix,
                          args.dirname,
                          args.pheno_file)


if __name__ == '__main__':
    main()
