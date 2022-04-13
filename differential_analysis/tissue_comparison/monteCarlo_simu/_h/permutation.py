"""
This script uses permutation testing to determine significance
of overlap between four brain regions.
"""
import pandas as pd
from random import sample
from functools import lru_cache

@lru_cache()
def get_mash():
    return pd.read_csv("../../_m/genes/lfsr_feature_4tissues.txt.gz",
                       sep='\t')


@lru_cache()
def get_sig_de(tissue):
    new_tissue = tissue.replace(" ", ".")
    return get_mash()[(get_mash()[new_tissue] < 0.05)].copy()


@lru_cache()
def get_regions():
    cc = get_sig_de("Caudate").Effect
    gg = get_sig_de("Dentate Gyrus").Effect
    dd = get_sig_de("DLPFC").Effect
    hh = get_sig_de("Hippocampus").Effect
    total = set(cc) | set(gg) | set(dd) | set(hh)
    return cc, gg, dd, hh, total


def perm_test(nc, ng, nd, nh, total, nmc):
    k = 0.0; diff = 112
    for i in range(nmc):
        cc = sample(total, nc)
        gg = sample(total, ng)
        dd = sample(total, nd)
        hh = sample(total, nh)
        k += diff <= len(set(cc) & set(gg) & set(dd) & set(hh))
    return k / nmc


def mc_perm_test(nmc):
    nc, ng, nd, nh, total = get_regions()
    return perm_test(len(nc), len(ng), len(nd), len(nh), sorted(total), nmc)


def main():
    perm = mc_perm_test(100000)
    print(f"Permutation p-value: {perm}")


if __name__ == "__main__":
    main()
