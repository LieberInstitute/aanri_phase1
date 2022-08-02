"""
This script prepares the eQTL results (tensorQTL) for mash modeling.
"""
import argparse
import pandas as pd

def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100,
                     usecols=["phenotype_id","variant_id","pval_g","pval_emt"])
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict(),
                       usecols=["phenotype_id","variant_id","pval_g","pval_emt"],
                       compression="gzip")


def extract_eqtls(feature):
    ## Load eQTLs for mashr
    ### Caudate
    print("Loading Caudate data!")
    cc_file = "../../../caudate/%s/tensoreqtl/_m/" % feature+\
        "LIBD_TOPMed_AA.cis_qtl_top_assoc.txt.gz"
    caudate = load_eqtl(cc_file)
    ### Dentate Gyrus
    print("Loading Dentate Gyrus data!")
    gg_file = "../../../dentateGyrus/%s/tensoreqtl/_m/" % feature+\
        "LIBD_TOPMed_AA.cis_qtl_top_assoc.txt.gz"
    gyrus = load_eqtl(gg_file)
    ### DLPFC
    print("Loading DLPFC data!")
    dd_file = "../../../dlpfc/%s/tensoreqtl/_m/" % feature+\
        "LIBD_TOPMed_AA.cis_qtl_top_assoc.txt.gz"
    dlpfc = load_eqtl(dd_file)
    ### Hippocampus
    print("Loading hippocampus data!")
    hh_file = "../../../hippocampus/%s/tensoreqtl/_m/" % feature+\
        "LIBD_TOPMed_AA.cis_qtl_top_assoc.txt.gz"
    hippo = load_eqtl(hh_file)
    return caudate, gyrus, dlpfc, hippo


def extract_dataframe(caudate, gyrus, dlpfc, hippo):
    ## Caudate
    df = pd.concat([caudate, gyrus, dlpfc, hippo])
    df["effect"] = df.phenotype_id + "_" + df.variant_id
    df.drop_duplicates(subset="effect")\
      .to_csv("top_eqtls_AA_interaction.txt.gz",
              sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    print("Load data!")
    cc, gg, dd, hh = extract_eqtls(args.feature)
    print("Extract strong set!")
    extract_dataframe(cc, gg, dd, hh)


if __name__=='__main__':
    main()
