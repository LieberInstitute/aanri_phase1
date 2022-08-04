"""
This script prepares the eQTL results (tensorQTL) for mash modeling.
"""
import argparse
import pandas as pd

def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100,
                     usecols=["phenotype_id","variant_id","slope","slope_se"])
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict(),
                       usecols=["phenotype_id","variant_id","slope","slope_se"],
                       compression="gzip")


def extract_eqtls(feature):
    ## Load eQTLs for mashr
    ### Caudate
    print("Loading Caudate data!")
    cc_file = "../../../caudate/%s/cis_analysis/_m/" % feature+\
        "LIBD_TOPMed_AA.nominal.txt.gz"
    caudate = load_eqtl(cc_file)
    ### Dentate Gyrus
    print("Loading Dentate Gyrus data!")
    gg_file = "../../../dentateGyrus/%s/cis_analysis/_m/" % feature+\
        "LIBD_TOPMed_AA.nominal.txt.gz"
    gyrus = load_eqtl(gg_file)
    ### DLPFC
    print("Loading DLPFC data!")
    dd_file = "../../../dlpfc/%s/cis_analysis/_m/" % feature+\
        "LIBD_TOPMed_AA.nominal.txt.gz"
    dlpfc = load_eqtl(dd_file)
    ### Hippocampus
    print("Loading hippocampus data!")
    hh_file = "../../../hippocampus/%s/cis_analysis/_m/" % feature+\
        "LIBD_TOPMed_AA.nominal.txt.gz"
    hippo = load_eqtl(hh_file)
    return caudate, gyrus, dlpfc, hippo


def extract_dataframe(caudate, gyrus, dlpfc, hippo, variable, label):
    ## Caudate
    dfc = caudate.loc[:, ["phenotype_id","variant_id",variable]]\
                 .rename(columns={variable: "Caudate"})
    ## Dentate Gyrus
    dfg = gyrus.loc[:, ["phenotype_id","variant_id",variable]]\
                 .rename(columns={variable: "Dentate_Gyrus"})
    # DLPFC
    dfd = dlpfc.loc[:, ["phenotype_id","variant_id",variable]]\
               .rename(columns={variable: "DLPFC"})
    ## Hippocampus
    dfh = hippo.loc[:, ["phenotype_id","variant_id",variable]]\
               .rename(columns={variable: "Hippocampus"})
    df = dfc.merge(dfg, on=["phenotype_id", "variant_id"])\
            .merge(dfd, on=["phenotype_id", "variant_id"])\
            .merge(dfh, on=["phenotype_id", "variant_id"])
    df.to_csv("%s_nominal_4tissues_AA.txt.gz" % label,
              sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    print("Load data!")
    cc, gg, dd, hh = extract_eqtls(args.feature)
    print("Subset for bhat!")
    extract_dataframe(cc, gg, dd, hh, "slope", "bhat")
    print("Subset for shat!")
    extract_dataframe(cc, gg, dd, hh, "slope_se", "shat")
    #print("Subset for pvalue!")
    #extract_dataframe(cc, gg, dd, hh, "pval_gi", "pval")


if __name__=='__main__':
    main()
