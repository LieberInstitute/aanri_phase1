"""
Examining the variation explained by ancestry-associated DE features.
"""

import numpy as np
import pandas as pd
import errno, os, argparse
from functools import lru_cache
from scipy.stats import linregress
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from rpy2.robjects import r, globalenv, pandas2ri

@lru_cache()
def get_ancestry():
    """
    Loads admixture ancestry from STRUCTURE.
    """
    fn = "../../../../../input/ancestry_structure/"+\
        "structure.out_ancestry_proportion_raceDemo_compare"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def get_pheno_data():
    """
    Load phenotype data.
    """
    fn = '../../../../../input/phenotypes/merged/_m/merged_phenotypes.csv'
    return pd.read_csv(fn, index_col=0)\
             .merge(get_ancestry(), left_on="BrNum", right_on="id")


@lru_cache()
def get_residualized(tissue, feature):
    '''
    Load residualization file.
    '''
    new_tissue = [tissue.lower() if tissue != "Dentate Gyrus" else "dentateGyrus"]
    fn = '../../../../%s/_m/%s/residualized_expression.tsv' % (new_tissue[0], feature)
    return pd.read_csv(fn, sep='\t', index_col=0).transpose()


@lru_cache()
def get_mash(feature):
    """
    Get local false sign rate (all)
    """
    fn = "../../_m/%s/lfsr_feature_4tissues.txt.gz" % feature
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_sig(tissue, feature, simu):
    """
    Get significant associations for specific brain region
    """
    df = get_mash(feature)[(get_mash(feature)[tissue.replace(" ", ".")]<0.05)]\
        .sort_values(tissue.replace(" ", "."))
    if simu == 0:
        return df
    else:
        return df.head(100)


@lru_cache()
def get_random(feature, simu):
    '''
    Select genes randomly.
    '''
    state = 13131313 + simu
    return get_mash(feature).sample(n=100, random_state=state)


@lru_cache()
def get_res_df(simu, feature, tissue, SIG=True):
    """
    Merge DE and residualized data after selecting # features.
    """
    if SIG:
        newList = list(get_sig(tissue, feature, simu).index)
    else:
        newList = list(get_random(feature, simu).index)
    return get_residualized(tissue, feature)[newList]


@lru_cache()
def get_deg_res_df(feature, tissue, num):
    """
    Merge DE and residualized data after selecting N features.
    """
    newList = list(get_sig(tissue, feature, 0)\
                   .sort_values(tissue.replace(" ", ".")).head(num).index)
    return get_residualized(tissue, feature)[newList]


def cal_pca(df):
    x = StandardScaler().fit_transform(df)
    pca = PCA(n_components=2).fit_transform(x)
    return pd.DataFrame(data=pca, columns=['PC1', 'PC2'], index=df.index)


def get_corr(dft):
    '''This calculates R^2 correlation via linear regression:
         - used to calculate relationship between 2 arrays
         - the arrays are principal components 1 or 2 (PC1, PC2) AND ancestry
         - calculated on a scale of 0 to 1 (with 0 being no correlation)
        Inputs:
            dft: Data frame with continuous variable and PCs
        Outputs:
          1. r2
          2. p-value, two-sided test
            - whose null hypothesis is that two sets of data are uncorrelated
          3. slope (beta): directory of correlations
    '''
    xx = dft.Eur; yy = dft.PC1
    slope, intercept, r_value, p_value, std_err = linregress(xx, yy)
    return r_value**2, p_value


def mkdir_p(directory):
    """
    Make a directory if it does not already exist.

    Input: Directory name
    """
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def get_pca_df(tissue, feature, num):
    '''
    new_pheno: This gets the correct size of samples using the the first two
               columns of residualized expression
      - the residualized expression data frame, has the correct samples
      - output new_pheno shape row numbers should be the same as res_df row numbers
    '''
    expr_res = get_res_df(num, feature, tissue, True)
    pheno_df = get_pheno_data()
    # Generate pheno data frame with correct samples
    new_pheno = pheno_df.merge(expr_res.iloc[:, 0:1], left_on="RNum",
                               right_index=True)\
                        .drop(expr_res.iloc[:, 0:1].columns, axis=1)\
                        .set_index("RNum")
    principalDf = cal_pca(expr_res)
    #get_explained_variance(tissue, expr_res)
    return pd.concat([principalDf, new_pheno], axis = 1)


def plotNsave(tissue, feature, fn, num):
    pandas2ri.activate()
    globalenv['df'] = get_pca_df(tissue, feature, num)
    globalenv['fn'] = fn
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.png', '.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = "Genetic ancestry (EA)"
    pp = ggpubr::ggscatter(df, x="Eur", y="PC1", add="reg.line", size=1,
                           xlab=xlab, panel.labs.font=list(face="bold"),
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.method="pearson", cor.coeff.args=list(label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15), ncol=4)
    save_plot(pp, fn, 6, 6)
    ''')


def pc_recursive_degs(tissue, feature):
    new_tissue = [tissue.lower() if tissue != "Dentate Gyrus" else "dentateGyrus"]
    pheno_df = get_pheno_data(); pvals = []; rsq = []; nums = []
    for num in range(2, 251):
        expr_res = get_deg_res_df(feature, tissue, num)
        new_pheno = pheno_df.merge(expr_res.iloc[:, 0:1],
                                   right_index=True, left_on="RNum")\
                            .drop(expr_res.iloc[:, 0:1].columns, axis=1)\
                            .set_index("RNum")
        dft = pd.concat([cal_pca(expr_res), new_pheno], axis = 1)
        r2, pval = get_corr(dft)
        nums.append(num); pvals.append(pval); rsq.append(r2)
    pd.DataFrame({"DEGs":nums, "PValue":pvals, "Rsq": rsq})\
        .sort_values("DEGs", ascending=True)\
        .to_csv("%s/variance_explained_rsq_250_de_%s.csv" %
                (new_tissue[0], feature), index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tissue', type=str)
    args=parser.parse_args()
    # Script calling
    tissue = args.tissue
    new_tissue = [tissue.lower() if tissue != "Dentate Gyrus" else "dentateGyrus"]
    mkdir_p(new_tissue[0])
    with open("%s/summarize_explained_variance.log" % new_tissue[0], mode="w") as f:
        print("%s, Explained Variance (%%) of PC1:" % tissue, file=f)
        for feature in ["genes", "transcripts", "exons", "junctions"]:
            r2, pval = get_corr(get_pca_df(tissue, feature, 0))
            print("All DEs\t%s\tr2: %0.2f%%\tp-value: %0.2e" %
                  (feature, r2 * 100, pval), file=f)
            plotNsave(tissue, feature, "%s/%s_de_pca_all" %
                      (new_tissue[0], feature), 0)
            r2, pval = get_corr(get_pca_df(tissue, feature, 100))
            print("Top 100 DEs\t%s\tr2: %0.2f%%\tp-value: %0.2e" %
                  (feature, r2 * 100, pval), file=f)
            plotNsave(tissue, feature, "%s/%s_de_pca_top100" %
                      (new_tissue[0], feature), 100)
            pc_recursive_degs(tissue, feature)


if __name__ == '__main__':
    main()
