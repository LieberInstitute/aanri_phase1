"""
Explore variation of ancestry PCs across multiple features.
"""

import numpy as np
import pandas as pd
import functools, argparse
from scipy.stats import linregress
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

@functools.lru_cache()
def get_deg(feature):
    '''
    Take DE genes obtained from limma-voom pipeline.
    '''
    fn = '../../_m/%s/diffExpr_EAvsAA_full.txt' % feature
    deg = pd.read_csv(fn, sep='\t', index_col=0).sort_values('adj.P.Val')
    return deg[(deg['adj.P.Val'] < 0.05)]


@functools.lru_cache()
def get_residualized(feature):
    '''
    Load residualization file.
    '''
    fn = '../../_m/%s/residualized_expression.tsv' % feature
    return pd.read_csv(fn, sep='\t', index_col=0).transpose()


@functools.lru_cache()
def get_ancestry():
    """
    Loads admixture ancestry from STRUCTURE.
    """
    fn = "../../../../input/ancestry_structure/"+\
        "structure.out_ancestry_proportion_raceDemo_compare"
    return pd.read_csv(fn, sep='\t')


@functools.lru_cache()
def get_pheno_data():
    """
    Load phenotype data.
    """
    fn = '../../../../input/phenotypes/_m/dg_phenotypes.csv'
    return pd.read_csv(fn, index_col=0)\
             .merge(get_ancestry(), left_on="BrNum", right_on="id")\
             .drop_duplicates(subset="BrNum")


@functools.lru_cache()
def get_deg_res_df(feature, num):
    """
    Merge DE and residualized data after selecting # features.
    """
    newList = list(get_deg(feature).sort_values("P.Value").head(num).index)
    return get_residualized(feature)[newList]


def cal_pca(df):
    x = StandardScaler().fit_transform(df)
    pca = PCA(n_components=2).fit_transform(x)
    return pd.DataFrame(data=pca, columns=['PC1', 'PC2'], index=df.index)


def get_pca_df(feature, num):
    '''
    new_pheno: This gets the correct size of samples using the the first two
               columns of residualized expression
      - the residualized expression data frame, has the correct samples
      - output new_pheno shape row numbers should be the same as res_df row numbers
    '''
    expr_res = get_deg_res_df(feature, num)
    pheno_df = get_pheno_data()
    # Generate pheno data frame with correct samples
    new_pheno = pheno_df.merge(expr_res.iloc[:, 0:1], right_index=True, left_index=True)\
                        .drop(expr_res.iloc[:, 0:1].columns, axis=1)
    principalDf = cal_pca(expr_res)
    get_explained_variance(expr_res)
    return pd.concat([principalDf, new_pheno], axis = 1)


def calculate_corr(xx, yy):
    '''This calculates R^2 correlation via linear regression:
         - used to calculate relationship between 2 arrays
         - the arrays are principal components 1 or 2 (PC1, PC2) AND ancestry
         - calculated on a scale of 0 to 1 (with 0 being no correlation)
        Inputs:
          x: array of variable of interest (continous or binary)
          y: array of PC
        Outputs:
          1. r2
          2. p-value, two-sided test
            - whose null hypothesis is that two sets of data are uncorrelated
          3. slope (beta): directory of correlations
    '''
    slope, intercept, r_value, p_value, std_err = linregress(xx, yy)
    return slope, r_value, p_value


def get_corr(dft):
    xx = dft.Eur
    yy = dft.PC1
    slope1, r_value1, p_value1 = calculate_corr(xx, yy)
    return r_value1**2, p_value1


def pc_recursive(feature):
    pheno_df = get_pheno_data()
    pvals = []; rsq = []; nums = []
    for num in range(2, 1001):
        expr_res = get_deg_res_df(feature, num)
        new_pheno = pheno_df.merge(expr_res.iloc[:, 0:1],
                                   right_index=True, left_on="RNum")\
                            .drop(expr_res.iloc[:, 0:1].columns, axis=1)\
                            .set_index("RNum")
        principalDf = cal_pca(expr_res)
        dft = pd.concat([principalDf, new_pheno], axis = 1)
        r2, pval = get_corr(dft)
        nums.append(num); pvals.append(pval); rsq.append(r2)
    pd.DataFrame({"DEGs":nums, "PValue":pvals, "Rsq": rsq})\
        .sort_values("Rsq", ascending=False)\
        .to_csv("%s_variance_explained_rsq_1000.csv" % feature, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    pc_recursive(args.feature)


if __name__ == '__main__':
    main()
