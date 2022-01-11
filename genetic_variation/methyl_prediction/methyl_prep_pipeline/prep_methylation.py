#!/usr/bin/env python

import functools
import pandas as pd
import os, argparse, re, errno

__author__ = 'KJ Benjamin'


@functools.lru_cache()
def get_phenotypes(pheno_file):
    """
    Get phenotypes for matching BrNums and expression.
    """
    return pd.read_csv(pheno_file, index_col=0)


@functools.lru_cache()
def get_expression(expr_file):
    """
    Get's the expression for all features. Recommend residualized expression.
    """
    return pd.read_csv(expr_file, sep='\t')


def mkdir_p(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def get_methylation(cpg_dir, gene_id):
    return pd.read_csv("%s/%s/cpg_methylation_data.csv" %
                       (cpg_dir, gene_id),
                       index_col=0)


def add_expression_n_save(expr_file, gene_file, dirname, pheno_file, cpg_dir):
    gene_id = re.search("ENSG\d+.\d+", gene_file).group(0)
    res_df = pd.merge(get_phenotypes(pheno_file),
                      get_expression(expr_file).loc[gene_id, :],
                      left_index=True, right_index=True)\
               .loc[:, ["BrNum", gene_id]]
    methyl_df = get_methylation(cpg_dir, gene_id).T
    if methyl_df.shape[1] != 0:
        mkdir_p(dirname)
        methyl_df.merge(res_df, left_index=True, right_on="BrNum")\
                 .drop_duplicates(gene_id)\
                 .drop(["BrNum"], axis=1)\
                 .to_csv("%s/methylation.csv" % dirname)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str)
    parser.add_argument('--cpg_dir', type=str)
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
                          args.cpg_dir)


if __name__ == '__main__':
    main()
