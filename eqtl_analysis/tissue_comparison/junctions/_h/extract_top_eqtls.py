"""
This script prepares the eQTL results (tensorQTL) for mash modeling.
"""
import argparse
import pandas as pd

def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100,
                     usecols=["phenotype_id","variant_id","pval_nominal"])
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict(),
                       usecols=["phenotype_id","variant_id","pval_nominal"],
                       compression="gzip")


def extract_eqtls(feature, fn):
    ## Load eQTLs for mashr
    ### Caudate
    print("Loading Caudate data!")
    cc_file = "../../../caudate/%s/cis_analysis/_m/%s" % (feature,fn)
    caudate = load_eqtl(cc_file)
    ### Dentate Gyrus
    print("Loading Dentate Gyrus data!")
    gg_file = "../../../dentateGyrus/%s/cis_analysis/_m/%s" % (feature,fn)
    gyrus = load_eqtl(gg_file)
    ### DLPFC
    print("Loading DLPFC data!")
    dd_file = "../../../dlpfc/%s/cis_analysis/_m/%s" % (feature,fn)
    dlpfc = load_eqtl(dd_file)
    ### Hippocampus
    print("Loading hippocampus data!")
    hh_file = "../../../hippocampus/%s/cis_analysis/_m/%s" % (feature,fn)
    hippo = load_eqtl(hh_file)
    return caudate, gyrus, dlpfc, hippo


def extract_dataframe(feature):
    fn1 = "LIBD_TOPMed_AA.signif_variants.txt.gz"
    fn2 = "LIBD_TOPMed_AA.conditional.txt.gz"
    print("Load data!")
    ## Top eQTL pairs
    cc1, gg1, dd1, hh1 = extract_eqtls(feature, fn1)
    # Conditional eQTL
    cc2, gg2, dd2, hh2 = extract_eqtls(feature, fn2)
    # Combine eQTL
    print("Extract strong set!")
    caudate = pd.concat([cc1, cc2]); gyrus = pd.concat([gg1, gg2])
    dlpfc = pd.concat([dd1, dd2]); hippo = pd.concat([hh1, hh2])
    # Combine brain regions
    df = pd.concat([caudate, gyrus, dlpfc, hippo])
    df["effect"] = df.phenotype_id + "_" + df.variant_id
    # Save files
    df.drop_duplicates(subset="effect")\
      .to_csv("top_eqtls_AA_nominal.txt.gz",sep='\t',index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    extract_dataframe(args.feature)


if __name__=='__main__':
    main()
