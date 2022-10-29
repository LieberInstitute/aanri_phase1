"""
This script runs tensorQTL in python.
"""

import argparse
import pandas as pd
from functools import lru_cache
from tensorqtl.post import get_significant_pairs

print(f"Pandas {pd.__version__}")

@lru_cache()
def get_permutation_results(prefix):
    return pd.read_csv("%s.genes.txt.gz" % (prefix), sep="\t", index_col=0)


@lru_cache()
def get_eqtl(prefix):
    return get_significant_pairs(get_permutation_results(prefix), prefix)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', type=str, default="AA")
    args=parser.parse_args()

    # Load data
    prefix = "LIBD_TOPMed_%s" % args.ancestry

    # Load and save eQTL results
    get_eqtl(prefix)\
        .to_csv("%s.signif_variants.txt.gz" % prefix, sep='\t', index=False)


if __name__ == "__main__":
    main()
