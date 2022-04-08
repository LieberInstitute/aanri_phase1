"""
Examining the variation explained by ancestry-associated DE isoform level.
"""

import numpy as np
import pandas as pd
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
    fn = "../../../../input/ancestry_structure/"+\
        "structure.out_ancestry_proportion_raceDemo_compare"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def get_pheno_data():
    """
    Load phenotype data.
    """
    fn = '../../../../input/phenotypes/merged/_m/merged_phenotypes.csv'
    return pd.read_csv(fn, index_col=0)\
             .merge(get_ancestry(), left_on="BrNum", right_on="id")


@lru_cache()
def get_dg_brnum():
    dg_brnum = get_pheno_data()[(get_pheno_data()["Region"] == "DentateGyrus")].BrNum.values
    return get_pheno_data()[(get_pheno_data()["BrNum"].isin(dg_brnum))].copy()


@lru_cache()
def get_residualized(tissue, feature):
    '''
    Load residualization file.
    '''
    new_tissue = [tissue.lower() if tissue != "Dentate Gyrus" else "dentateGyrus"]
    fn = '../../../%s/_m/%s/residualized_expression.tsv' % (new_tissue[0], feature)
    return pd.read_csv(fn, sep='\t', index_col=0).transpose()


@lru_cache()
def get_mash(feature):
    """
    Get local false sign rate (all)
    """
    fn = "../../_m/%s/lfsr_feature_4tissues.txt.gz" % feature
    return pd.read_csv(fn, sep='\t', index_col=0)


@lru_cache()
def get_sig(tissue, feature):
    """
    Get significant associations for specific brain region
    """
    df = get_mash(feature)[(get_mash(feature)[tissue.replace(" ", ".")]<0.05)]\
        .sort_values(tissue.replace(" ", "."))
    return df


@lru_cache()
def get_res_df(tissue, feature):
    """
    Merge DE and residualized data after selecting # features.
    """
    return get_residualized(tissue,feature)[list(get_sig(tissue,feature).index)]


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
    '''
    xx = dft.Eur; yy = dft.PC1
    slope, intercept, r_value, p_value, std_err = linregress(xx, yy)
    return r_value**2, p_value


def get_pca_df(tissue, feature):
    '''
    new_pheno: This gets the correct size of samples using the the first two
               columns of residualized expression
      - the residualized expression data frame, has the correct samples
      - output new_pheno shape row numbers should be the same as res_df row numbers
    '''
    expr_res = get_res_df(tissue, feature)
    pheno_df = get_dg_brnum()
    # Generate pheno data frame with correct samples
    new_df = pheno_df.merge(expr_res, left_on="RNum", right_index=True)
    new_pheno = new_df.drop(expr_res.columns, axis=1).set_index("RNum")
    new_expr = new_df.drop(pheno_df.drop("RNum", axis=1).columns, axis=1)\
                     .set_index("RNum")
    principalDf = cal_pca(new_expr)
    #get_explained_variance(tissue, pd.concat([principalDf, new_pheno], axis=1))
    return pd.concat([principalDf, new_pheno], axis = 1)


def plotNsave(tissue, feature, fn):
    pandas2ri.activate()
    globalenv['df'] = get_pca_df(tissue, feature)
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


def main():
    with open("summarize_explained_variance.log", mode="w") as f:
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            new_tissue = [tissue.lower() if tissue != "Dentate Gyrus" else "dentateGyrus"]
            print("%s, Explained Variance (%%) of PC1:" % tissue, file=f)
            for feature in ["genes", "transcripts", "exons", "junctions"]:
                r2, pval = get_corr(get_pca_df(tissue, feature))
                plotNsave(tissue, feature, "%s_deg_pca_%s" %
                          (new_tissue[0], feature))
                print("%s\tr2: %0.2f%%\tp-value: %0.2e" %
                      (feature, r2 * 100, pval), file=f)


if __name__ == '__main__':
    main()
