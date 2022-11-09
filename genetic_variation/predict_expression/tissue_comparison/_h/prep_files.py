"""
This script prepares DE analysis for mash modeling.
"""
import pandas as pd
import argparse, session_info
from functools import lru_cache

def load_de(filename):
    return pd.read_csv(filename, sep='\t', index_col=0)


@lru_cache()
def extract_de(feature):
    ## Load DEs for mashr
    ### Caudate
    cc_file = "../../caudate/_m/%s/diffExpr_EAvsAA_full.txt" % feature
    caudate = load_de(cc_file)
    ### Dentate Gyrus
    gyrus_file = "../../dentateGyrus/_m/%s/diffExpr_EAvsAA_full.txt" %feature
    gyrus = load_de(gyrus_file)
    ### DLPFC
    dlpfc_file = "../../dlpfc/_m/%s/diffExpr_EAvsAA_full.txt" % feature
    dlpfc = load_de(dlpfc_file)
    ### Hippocampus
    hippo_file = "../../hippocampus/_m/%s/diffExpr_EAvsAA_full.txt" % feature
    hippo = load_de(hippo_file)
    return caudate, gyrus, dlpfc, hippo


def extract_dataframe(variable, label, feature):
    caudate, gyrus, dlpfc, hippo = extract_de(feature)
    ## Caudate
    df1 = caudate.loc[:, [variable]].rename(columns={variable: "Caudate"})
    ## Dentate Gyrus
    df2 = gyrus.loc[:, [variable]].rename(columns={variable: "Dentate Gyrus"})
    ## DLPFC
    df3 = dlpfc.loc[:, [variable]].rename(columns={variable: "DLPFC"})
    ## Hippocampus
    df4 = hippo.loc[:, [variable]].rename(columns={variable: "Hippocampus"})
    df = df1.merge(df2, left_index=True, right_index=True)\
            .merge(df3, left_index=True, right_index=True)\
            .merge(df4, left_index=True, right_index=True)\
            .reset_index().rename(columns={"index": "Feature"})
    df.to_csv("%s/%s_de_4tissues.tsv" % (feature, label), sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    extract_dataframe("logFC", "bhat", args.feature)
    extract_dataframe("SE", "shat", args.feature)
    ## Session information
    session_info.show()


if __name__=='__main__':
    main()
