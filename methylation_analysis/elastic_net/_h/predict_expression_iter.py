"""
This script runs tensorQTL in python.
"""
import pandas as pd
from pyhere import  here
from functools import lru_cache
import argparse, errno, os, session_info
from tensorqtl import read_phenotype_bed
from sklearn.model_selection import KFold
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import explained_variance_score
from sklearn.metrics import r2_score, mean_squared_error


def map_tissue(tissue):
    return {"caudate": "Caudate", "dlpfc": "DLPFC",
            "dentateGyrus": "Dentate_Gyrus", "hippocampus": "Hippocampus"}[tissue]


def map_feature(feature):
    return {"genes": "Gene", "transcripts": "Transcript",
            "exons": "Exon", "junctions": "Junction"}[feature]


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


@lru_cache()
def load_de(tissue, feature):
    de_file = here("differential_analysis/tissue_comparison/summary_table/_m/",
                   "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    de_df = pd.read_csv(de_file, sep='\t')
    return de_df[(de_df["Tissue"] == map_tissue(tissue).replace("_", " ")) &
                 (de_df["Type"] == map_feature(feature))].copy()


def get_phenotype(tissue, feature):
    expr_bed = here("eqtl_analysis/%s/%s/"% (tissue, feature),
                    "normalize_expression/_m/",
                    "%s.expression.bed.gz" % feature)
    return read_phenotype_bed(str(expr_bed))


@lru_cache()
def subset_phenotypes(tissue, feature):
    phenotype_df, phenotype_pos_df = get_phenotype(tissue, feature)
    shared_features = set(load_de(tissue, feature).Effect) & \
        set(phenotype_df.index)
    phenotype_pos_df["gene_id"] = phenotype_pos_df.index
    phenotype_df = phenotype_df.loc[shared_features].sort_index()
    phenotype_pos_df = phenotype_pos_df.loc[shared_features].sort_index()
    return phenotype_df, phenotype_pos_df


@lru_cache()
def get_methylation(chrom):
    methyl_df = pd.read_csv("../_m/methylation/methylation_wgbs_%s.tsv" % chrom,
                            sep='\t', index_col=0)
    annot_df = pd.read_csv("../_m/methylation/annotation_CpGs_wgbs_%s.tsv" % chrom,
                           sep='\t', index_col=0)
    annot_df = annot_df.loc[:, ["seqnames", "start"]]\
                       .rename(columns={"seqnames":"chrom",
                                        "start":"pos"})
    return methyl_df, annot_df


def load_data(tissue, feature, ID, window=500000):
    phenotype_df, phenotype_pos_df = subset_phenotypes(tissue, feature)
    (chrom, tss, geneid) = phenotype_pos_df.iloc[ID, :]
    methyl_df, annot_df = get_methylation(chrom)
    mask = (annot_df["pos"] <= (tss+window)) & (annot_df["pos"] >= (tss-window))
    annot_df = annot_df[(annot_df["chrom"] == chrom) & mask].copy()
    cpgs_ids = annot_df.index
    samples = set(methyl_df.columns) & set(phenotype_df.columns)
    return methyl_df.loc[cpgs_ids,samples].T, phenotype_df.loc[geneid,samples]


def model_preformance(estimator, x_train, x_test, y_train, y_test, fold_num):
    estimator.fit(x_train, y_train)
    labels_pred = estimator.predict(x_test)
    pred_df = pd.DataFrame(labels_pred, index=x_test.index,
                           columns=[y_test.name])
    weights_df = pd.DataFrame(estimator.coef_,
                              index=x_test.columns,
                              columns=["fold_%d" % fold_num])
    output = dict()
    output['feature'] = y_test.name
    output['n_features'] = len(x_train.columns)
    output['test_score_r2'] = r2_score(y_test, labels_pred)
    output['test_score_mse'] = mean_squared_error(y_test, labels_pred)
    output['test_score_evar'] = explained_variance_score(y_test, labels_pred,
                                                         multioutput='uniform_average')
    return output, pred_df, weights_df


def compute_expression(df, X, Y):
    return pd.DataFrame({Y.name:X.dot(df.mean(axis=1).values)},
                        index=X.index)


def predict_expression(tissue, feature, sge_id, nsplit):
    """
    Function to run each feature. This is for parallelization.
    """
    regr = ElasticNetCV(cv=5, random_state=20221003)
    kf = KFold(n_splits=nsplit, shuffle=True, random_state=13)
    fields = ['feature', 'n_features', 'test_score_r2',
              'test_score_mse', 'test_score_evar']
    outdir = "%s" % (feature); mkdir_p(outdir)
    X, Y = load_data(tissue, feature, sge_id, window=500000)
    if X.shape[1] != 0:
        fold = 0; df_dict = pd.DataFrame(); weights_df = pd.DataFrame()
        with open('%s/enet_test_metrics_%d.txt' % (outdir, sge_id), 'w') as f:
            print("\t".join(['fold'] + fields), file=f, flush=True)        
            for train_index, test_index in kf.split(X):
                X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
                Y_train, Y_test = Y[train_index], Y[test_index]
                o, _, weights = model_preformance(regr, X_train, X_test,
                                                  Y_train, Y_test, fold)
                # df_dict = pd.concat([df_dict, pred_df], axis=0)
                weights_df = pd.concat([weights_df, weights], axis=1)
                print("\t".join([str(fold)] + [str(o[x]) for x in fields]),
                      flush=True, file=f)
                fold += 1
            # df_dict.to_csv("%s/predicted_expression_%d.txt" % (outdir, sge_id),
            #                sep='\t', index=True)
        compute_expression(weights_df, X, Y)\
            .to_csv("%s/predicted_expression_%d.txt" % (outdir, sge_id),
                    sep='\t', index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tissue', type=str)
    parser.add_argument('--feature', type=str)
    parser.add_argument('--nsplit', type=int)
    parser.add_argument('--start', type=int)
    parser.add_argument('--end', type=int)
    args=parser.parse_args()

    # Predict expression with elastic net
    tissue = args.tissue; feature = args.feature; nsplit = args.nsplit
    start = args.start; end = args.end
    pred_it = ((tissue, feature, k, nsplit) for k in range(start, end))
    for tt,ff,ii,nn in pred_it:
        predict_expression(tt, ff, ii, nn)
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
