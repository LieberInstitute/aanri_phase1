## This script calculates exact p-values for
## Spearman correlation.

import session_info
import polars as pl
from scipy.stats import spearmanr

def load_dmrs():
    fn = "../_m/dmr_global_local_combined_3regions.csv"
    return pl.scan_csv(fn)


def select_region(region):
    return load_dmrs()\
        .filter(pl.col("region") == region)


def get_sig_dmrs(region):
    return select_region(region)\
        .filter((pl.col("fdr_global") < 0.05) |
                (pl.col("fdr_local") < 0.05))


def corr_beta(region, SIG=False):
    if SIG:
        df = get_sig_dmrs(region).collect()
    else:
        df = select_region(region).collect()
    return spearmanr(df.select(pl.col("beta_global")),
                     df.select(pl.col("beta_local")))


def main():
    ## Correlated effect sizes
    with open("rho_statistics.log", "w") as f:
        for SIG in [False, True]:
            for region in ["Caudate", "DLPFC", "Hippocampus"]:
                rho, pval = corr_beta(region, SIG)
                print(f"Significant ({SIG})\t",
                      f"{region}:\t rho > {rho: .3f}\t",
                      f"p-value < {pval: .2e}",
                      file=f)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
