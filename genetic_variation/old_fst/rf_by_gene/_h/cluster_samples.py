import functools
import pandas as pd


@functools.lru_cache()
def get_fam():
    fam_file = "/ceph/projects/brainseq/genotype/download/topmed/"+\
        "convert2plink/filter_maf_01/_m/LIBD_Brain_TopMed.fam"
    return pd.read_csv(fam_file, sep='\t', header=None,
                       names=["FID", "IID", "V2", "V3", "Sex", "Pheno"])

@functools.lru_cache()
def generate_cluster():
    ancestry_file = "../../../../input/ancestry_structure/structure."+\
        "out_ancestry_proportion_raceDemo_compare"
    df = get_fam().merge(pd.read_csv(ancestry_file, sep='\t'),
                         left_on="FID", right_on="id")
    return df.loc[:, ["FID", "IID", "group"]].drop_duplicates(subset="FID")


@functools.lru_cache()
def get_phenotypes():
    pheno_file = "../../../../input/phenotypes/merged/_m/individual_phenotypes.csv"
    return pd.read_csv(pheno_file)


@functools.lru_cache()
def get_expr_data():
    dt = generate_cluster().merge(get_phenotypes(), left_on="FID",
                                  right_on="BrNum")
    return dt.loc[(dt["Age"] > 17) & (dt["Dx"] == "Control"),
                  ["FID", "IID", "group"]]


def main():
    ## Cluster names
    get_expr_data()\
        .to_csv("within_clusters.txt", sep='\t', index=False)
    ## Samples to include
    get_expr_data()\
        .loc[:, ["FID", "IID"]]\
        .to_csv("samples_keep.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
