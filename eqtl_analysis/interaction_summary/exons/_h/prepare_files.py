"""
This script prepares the eQTL results (tensorQTL) for mash modeling.
"""
import argparse
import pandas as pd

def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100,
                     usecols=["phenotype_id","variant_id","b_gi","b_gi_se"])
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict(),
                       usecols=["phenotype_id","variant_id","b_gi","b_gi_se"],
                       compression="gzip")


def extract_eqtls(feature):
    ## Load eQTLs for mashr
    ### Caudate
    cc_file = "../../../caudate/%s/interaction_model/_m/" % feature+\
        "LIBD_TOPMed_AA.interaction.txt.gz"
    caudate = load_eqtl(cc_file)
    ### Dentate Gyrus
    gg_file = "../../../dentateGyrus/%s/interaction_model/_m/" % feature+\
        "LIBD_TOPMed_AA.interaction.txt.gz"
    gyrus = load_eqtl(gg_file)
    ### DLPFC
    dd_file = "../../../dlpfc/%s/interaction_model/_m/" % feature+\
        "LIBD_TOPMed_AA.interaction.txt.gz"
    dlpfc = load_eqtl(dd_file)
    ### Hippocampus
    hh_file = "../../../hippocampus/%s/interaction_model/_m/" % feature+\
        "LIBD_TOPMed_AA.interaction.txt.gz"
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
    df.to_csv("%s_interaction_4tissues_AA.txt.gz" % label,
              sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    cc, gg, dd, hh = extract_eqtls(args.feature)
    extract_dataframe(cc, gg, dd, hh, "b_gi", "bhat")
    extract_dataframe(cc, gg, dd, hh, "b_gi_se", "shat")


if __name__=='__main__':
    main()
