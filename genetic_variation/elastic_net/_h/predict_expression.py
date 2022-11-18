"""
This script runs tensorQTL in python.
"""
import pandas as pd
from pyhere import  here
from functools import lru_cache
from sklearn.model_selection import KFold
from sklearn.linear_model import ElasticNetCV
import argparse, torch, errno, os, session_info
from sklearn.metrics import explained_variance_score
from tensorqtl import read_phenotype_bed, genotypeio
from tensorqtl.core import impute_mean, calculate_maf
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
    return phenotype_df.loc[shared_features], phenotype_pos_df.loc[shared_features]


@lru_cache()
def get_genotypes(tissue, feature):
    plink_prefix_path = here("eqtl_analysis/%s/%s/" % (tissue, feature),
                             "plink_format/_m/genotypes")
    pr = genotypeio.PlinkReader(str(plink_prefix_path))
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(fam_id="fid"), variant_df


def load_data(tissue, feature, ID, window=500000, maf_threshold=0.01):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    genotype_df, variant_df = get_genotypes(tissue, feature)
    phenotype_df, phenotype_pos_df = subset_phenotypes(tissue, feature)
    (geneid, chrom, tss) = phenotype_pos_df.reset_index().iloc[ID]
    mask = (variant_df["pos"] <= (tss+window)) & (variant_df["pos"] >= (tss-window))
    variant_df = variant_df[(variant_df["chrom"] == chrom) & mask].copy()
    variant_ids = variant_df.index
    genotypes = genotype_df.loc[variant_ids,:]
    genotypes_t = torch.tensor(genotypes.values, dtype=torch.float).to(device)
    impute_mean(genotypes_t)
    maf_t = calculate_maf(genotypes_t)
    mask_t = maf_t >= maf_threshold
    genotypes_t = genotypes_t[mask_t]
    mask = mask_t.cpu().numpy().astype(bool)
    variant_ids = variant_ids[mask]
    return genotype_df.loc[variant_ids,:].T, phenotype_df.loc[geneid,:]


def model_preformance(estimator, x_train, x_test, y_train, y_test):
    estimator.fit(x_train, y_train)
    labels_pred = estimator.predict(x_test)
    pred_df = pd.DataFrame(labels_pred, index=x_test.index,
                           columns=[y_test.name])
    output = dict()
    output['feature'] = y_test.name
    output['n_features'] = len(x_train.columns)
    output['test_score_r2'] = r2_score(y_test, labels_pred)
    output['test_score_mse'] = mean_squared_error(y_test, labels_pred)
    output['test_score_evar'] = explained_variance_score(y_test, labels_pred,
                                                         multioutput='uniform_average')
    return output, pred_df


def predict_expression(tissue, feature, sge_id, nsplit):
    """
    Function to run each feature. This is for parallelization.
    """
    regr = ElasticNetCV(cv=5, random_state=20221003)
    kf = KFold(n_splits=nsplit, shuffle=True, random_state=13)
    fields = ['feature', 'n_features', 'test_score_r2',
              'test_score_mse', 'test_score_evar']
    outdir = "%s/%s" % (tissue, feature); mkdir_p(outdir)
    X, Y = load_data(tissue, feature, sge_id, window=500000, maf_threshold=0.01)
    if X.shape[1] != 0:
        fold = 0; df_dict = pd.DataFrame()        
        with open('%s/enet_test_metrics_%d.txt' % (outdir, sge_id), 'w') as f:
            print("\t".join(['fold'] + fields), file=f, flush=True)        
            for train_index, test_index in kf.split(X):
                X_train, X_test = X.iloc[train_index, :], X.iloc[test_index, :]
                Y_train, Y_test = Y[train_index], Y[test_index]
                o, pred_df = model_preformance(regr, X_train, X_test,
                                               Y_train, Y_test)
                df_dict = pd.concat([df_dict, pred_df], axis=0)
                print("\t".join([str(fold)] + [str(o[x]) for x in fields]),
                      flush=True, file=f)
            fold += 1
        df_dict.to_csv("%s/predicted_expression_%d.txt" % (outdir, sge_id),
                       sep='\t', index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tissue', type=str)
    parser.add_argument('--feature', type=str)
    parser.add_argument('--sge_id', type=int)
    parser.add_argument('--nsplit', type=int)
    args=parser.parse_args()

    # Predict expression with elastic net
    tissue = args.tissue; feature = args.feature
    sge_id = args.sge_id; nsplit = args.nsplit
    predict_expression(tissue, feature, sge_id-1, nsplit)
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
