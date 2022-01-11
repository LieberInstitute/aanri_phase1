#!/usr/bin/env python

import functools
import pandas as pd
import os, argparse, re, errno

__author__ = 'KJ Benjamin'


@functools.lru_cache()
def get_phenotypes(pheno_file):
    """
    Get phenotypes for matching BrNums and expression.
    ------
    Input: Phenotype file
    """
    return pd.read_csv(pheno_file, index_col=0)


@functools.lru_cache()
def get_expression(expr_file):
    """
    Get's the expression for all features. Recommend residualized expression.
    ------
    Input: residualized expression file
    """
    return pd.read_csv(expr_file, sep='\t')


def mkdir_p(directory):
    """
    Make new directory if it does not exist. This is similar to mkdir -p.
    -----
    Input: directory name
    """
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def get_methylation(cpg_dir, gene_id):
    """
    Get CpG data that has already been prepared with onehot encoding.
    -------
    Inputs:
    _______

    cpg_dir: CpG prepared data directory PATH
    gene_id: Gene ID for analysis
    """
    return pd.read_csv("%s/%s/cpg_methylation_data.csv" %
                       (cpg_dir, gene_id),
                       index_col=0)


def get_snps(snp_dir, gene_id):
    """
    Get SNP data that has already been prepared with onehot encoding.
    -------
    Inputs:
    _______

    snp_dir: SNP onehot encoding directory PATH
    gene_id: Gene ID for analysis
    """
    gene = gene_id.replace(".", "_")
    return pd.read_csv("%s/%s/snps_onehot.csv" % (snp_dir, gene),
                       index_col=0).drop(gene_id, axis=1)


def add_expression_n_save(expr_file, gene_file, dirname, pheno_file,
                          cpg_dir, snp_dir):
    gene_id = re.search("ENSG\d+.\d+", gene_file).group(0)
    res_df = pd.merge(get_phenotypes(pheno_file),
                      get_expression(expr_file).loc[gene_id, :],
                      left_index=True, right_index=True)\
               .loc[:, ["BrNum", gene_id]]
    try:
        methyl_df = get_methylation(cpg_dir, gene_id)
        snp_df = get_snps(snp_dir, gene_id)
    except FileNotFoundError:
        methyl_df = pd.DataFrame()
    ## Combine and save
    if methyl_df.shape[1] != 0:
        mkdir_p(dirname)
        snp_df.merge(methyl_df, left_index=True, right_index=True)\
              .merge(res_df, left_index=True, right_index=True)\
              .drop(["BrNum"], axis=1)\
              .to_csv("%s/snps_methylation.csv" % dirname)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str)
    parser.add_argument('--cpg_dir', type=str)
    parser.add_argument('--snp_dir', type=str)
    parser.add_argument('--expr_file', type=str)
    parser.add_argument('--gene_file', type=str)
    parser.add_argument('--pheno_file', type=str)
    args=parser.parse_args()
    ## Main section
    os.environ['NUMEXPR_MAX_THREADS'] = '1'
    add_expression_n_save(args.expr_file,
                          args.gene_file,
                          args.dirname,
                          args.pheno_file,
                          args.cpg_dir,
                          args.snp_dir)


if __name__ == '__main__':
    main()
