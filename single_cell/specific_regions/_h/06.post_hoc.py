#### This script performs pairwise comparisons of means
#### for glial cell compositions

import session_info
import pandas as pd
from scipy.stats import tukey_hsd
from statsmodels.stats.multitest import fdrcorrection

def load_data(dataset):
    fn = f"transformed_proportions.{dataset.lower()}.tsv"
    return pd.read_csv(fn, sep='\t')


def subset_data(dataset):
    df = load_data(dataset)
    dlpfc = df[(df["region"] == "DLPFC")].copy()
    hpc = df[(df["region"] == "HPC")].copy()
    nac = df[(df["region"] == "NAC")].copy()
    return dlpfc, hpc, nac


def cal_tukey(dataset):
    dlpfc, hpc, nac = subset_data(dataset)
    clusters = dlpfc.clusters.unique()
    df = pd.DataFrame()
    for clust in clusters:
        res = tukey_hsd(dlpfc[(dlpfc["clusters"] == clust)].Freq,
                        hpc[(hpc["clusters"] == clust)].Freq,
                        nac[(nac["clusters"] == clust)].Freq)
        stats = pd.DataFrame(res.statistic, index=["DLPFC", "HPC", "NAc"])\
                  .rename(columns={0:"DLPFC", 1:"HPC", 2:"NAc"})
        pvals = pd.DataFrame(res.pvalue, index=["DLPFC", "HPC", "NAc"])\
                  .rename(columns={0:"DLPFC", 1:"HPC", 2:"NAc"})
        dx1 = pd.melt(pvals.reset_index(), id_vars=["index"],
                      var_name="region", value_name="pvalue")
        dx2 = pd.melt(stats.reset_index(), id_vars=["index"],
                      var_name="region", value_name="statistic")
        dx = pd.merge(dx1, dx2, on=["index", "region"])\
               .rename(columns={"index":"region_1", "region":"region_2"})
        dx = dx[~(dx["region_1"] == dx["region_2"])]\
            .sort_values("region_1")\
            .drop_duplicates(subset="pvalue")
        dx["cluster"] = clust
        df = pd.concat([df, dx], axis=0)
    _, fdr = fdrcorrection(df.pvalue)
    df["fdr"] = fdr
    return df


def main():
    ## Tukey HSD
    for dataset in ["Astro", "Oligo", "Micro"]:
        fn = f"propeller_proportion.{dataset.lower()}.tukey_hsd.tsv"
        cal_tukey(dataset).to_csv(fn, sep='\t', index=False)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
