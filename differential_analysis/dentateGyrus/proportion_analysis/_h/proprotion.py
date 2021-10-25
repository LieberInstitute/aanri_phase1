"""
This script explores the proportion of DE for all features.
"""
import functools
import pandas as pd

@functools.lru_cache()
def get_de(feature):
    fn = "../../_m/%s/diffExpr_EAvsAA_full.txt" % feature
    return pd.read_csv(fn, sep='\t', index_col=0)


@functools.lru_cache()
def get_sig_de(feature):
    return get_de(feature)[(get_de(feature)["adj.P.Val"] < 0.05)]


@functools.lru_cache()
def get_prop(feature):
    return get_sig_de(feature).shape[0] / get_de(feature).shape[0]


def main():
    out = "Proprotion of Ancestry DE:\n"+\
        "Gene - {:.1%}\n".format(get_prop("genes"))+\
        "Tx   - {:.1%}\n".format(get_prop("transcripts"))+\
        "Exon - {:.1%}\n".format(get_prop("exons"))+\
        "Jxn  - {:.1%}".format(get_prop("junctions"))
    with open("summary.log", "w") as f:
        print(out, file=f)


if __name__ == '__main__':
    main()
