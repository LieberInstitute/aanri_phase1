#!/usr/bin/python
"""
This was developed on Python 3.7+.
"""

import dRFEtools
import numpy as np
import pandas as pd
from matplotlib import rcParams
import errno, os, functools, argparse, re
from sklearn.model_selection import KFold
from rpy2.robjects import r, pandas2ri, globalenv
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import explained_variance_score
from sklearn.metrics import r2_score, mean_squared_error

@functools.lru_cache()
def get_phenotypes(pheno_file):
    """
    Get the phenotypes to use for stratified k-fold.
    """
    return pd.read_csv(pheno_file, index_col=0)


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


def run_enet(X_train, X_test, Y_train, Y_test, fold_num, outdir):
    """
    This function applys random forest with dynamic recusive
    feature elimination. Ranking features is default.

    Input: Training and test sets, fold number, and output directory name.

    Output:
    Text file: Ranked list of features
    Text file: Test prediction for model
    DataFrame: Data frame with feature elimination information
    Dictionary: Results of dynamic recursive feature elimination
    """
    # Apply elastic net with alpha tuning (default in dRFEtools)
    regr = dRFEtools.ElasticNet_RFE(n_jobs=2)
    features = X_train.columns
    d, pfirst = dRFEtools.dev_rfe(regr, X_train.values, Y_train.values,
                                  features, fold_num, out_dir=outdir,
                                  elimination_rate=0.2)
    df_elim = pd.DataFrame([{'fold':fold_num,
                             'n features':k,
                             'R2 Score':d[k][1],
                             'Mean Square Error':d[k][2],
                             'Explained Variance':d[k][3]} for k in d.keys()])
    # dRFEtools.plot_r2(d, fold_num, outdir)
    # dRFEtools.plot_mse(d, fold_num, outdir)
    # dRFEtools.plot_evar(d, fold_num, outdir)
    n_features_max = max(d, key=lambda x: d[x][1])
    # Max features based on lowess curve
    try:
        # Max features based on lowess curve
        n_features,_ = dRFEtools.extract_max_lowess(d, multi=False)
        # Redundant features based on lowess curve
        n_redundant,_ = dRFEtools.extract_redundant_lowess(d, multi=False)
        dRFEtools.plot_with_lowess_vline(d, fold_num, outdir,
                                         classify=False, multi=False)
    except ValueError:
        n_features = n_features_max
        n_redundant = n_features
    # Fit model
    regr.fit(X_train, Y_train)
    all_fts = regr.predict(X_test)
    regr.fit(X_train.values[:, d[n_redundant][4]], Y_train)
    labels_pred_redundant = regr.predict(X_test.values[:, d[n_redundant][4]])
    regr.fit(X_train.values[:, d[n_features][4]], Y_train)
    labels_pred = regr.predict(X_test.values[:, d[n_features][4]])
    # Output test predictions
    pd.DataFrame({'fold': fold_num,
                  'real': Y_test,
                  'predict_all': all_fts,
                  'predict_max': labels_pred,
                  'predict_redundant': labels_pred_redundant})\
      .to_csv('%s/test_predictions.txt' % outdir,
              sep='\t', mode='a', index=True,
              header=True if fold_num == 0 else False)
    # Save output data
    output = dict()
    output['n_features'] = n_features
    output['n_redundant'] = n_redundant
    output['n_features_all_features'] = pfirst[0]
    output['train_dev_score_r2_all_features'] = pfirst[1]
    output['train_dev_score_mse_all_features'] = pfirst[2]
    output['train_dev_score_evar_all_features'] = pfirst[3]
    output['test_score_r2'] = r2_score(Y_test, labels_pred)
    output['test_score_mse'] = mean_squared_error(Y_test, labels_pred)
    output['test_score_evar'] = explained_variance_score(Y_test, labels_pred,
                                                         multioutput='uniform_average')
    output['test_score_r2_redundant'] = r2_score(Y_test, labels_pred_redundant)
    output['test_score_mse_redundant'] = mean_squared_error(Y_test, labels_pred_redundant)
    output['test_score_evar_redundant'] = explained_variance_score(Y_test, labels_pred_redundant,
                                                                   multioutput='uniform_average')
    return output, df_elim


def run_by_feature(args):
    """
    Function to run each feature. This is for parallelization.
    """
    mkdir_p(args.tissue)
    dirname = re.search("(ENSG\d+_\d+)", args.fn).group(0)
    outdir = '%s/%s' % (args.tissue, dirname)
    gene_id = dirname.replace("_", ".")
    mkdir_p(outdir)
    df = pd.read_csv("%s/snps_onehot.csv" % args.fn, index_col=0)
    X = df.drop([gene_id], axis=1)
    Y = df.loc[:, gene_id]
    y = pd.DataFrame(Y).merge(get_phenotypes(), left_index=True,
                              right_index=True).Race.astype('category').cat.codes
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=13)
    skf.get_n_splits(X, y)
    fields = ['n_features_all_features', 'train_dev_score_mse_all_features',
              'train_dev_score_r2_all_features', 'train_dev_score_evar_all_features',
              'n_features', 'test_score_r2', 'test_score_mse',
              'test_score_evar', 'n_redundant', 'test_score_r2_redundant',
              'test_score_mse_redundant', 'test_score_evar_redundant']
    fold = 0
    df_dict = pd.DataFrame()
    with open('%s/enet_rfe_10folds.txt' % (outdir), 'w') as f:
        print("\t".join(['fold'] + fields), file=f, flush=True)
        for train_index, test_index in skf.split(X, y):
            #print("TRAIN:", train_index, "\nTEST:", test_index, "\nFold:", fold)
            X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
            Y_train, Y_test = Y[train_index], Y[test_index]
            o, df_elim = run_enet(X_train, X_test, Y_train, Y_test, fold, outdir)
            df_dict = pd.concat([df_dict, df_elim], axis=0)
            print("\t".join([str(fold)] + [str(o[x]) for x in fields]),
                  flush=True, file=f)
            fold += 1
        df_dict.to_csv("%s/feature_elimination_allFolds_metrics.txt" % outdir,
                       sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pheno_file', type=str)
    parser.add_argument('--tissue', type=str)
    parser.add_argument('--fn', type=str)
    args=parser.parse_args()
    os.environ['NUMEXPR_MAX_THREADS'] = '1'
    rcParams.update({'figure.max_open_warning': 0})
    run_by_feature(args)


if __name__ == '__main__':
    main()
