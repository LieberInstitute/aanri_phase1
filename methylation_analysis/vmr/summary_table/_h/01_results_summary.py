## This script examines the average metrics for DMRs
## (global and local) ancestry-associated.

import session_info
import polars as pl

def load_dmrs(tissue, LOCAL):
    if LOCAL:
        fn = f"../../_m/aaOnly/{tissue.lower()}_dmr_local.csv"
        ancestry = "local"
        df = pl.scan_csv(fn, sep="\t")
    else:
        fn = f"../../_m/aaOnly/{tissue.lower()}_dmr.csv"
        ancestry = "global"
        df = pl.scan_csv(fn, sep=",")
    return df.select([
               pl.col("*"),
               pl.when(pl.col("beta") > 0).then("hyper").otherwise("hypo").alias("dmr_status")])\
             .filter(pl.col("fdr") < 0.05)\
             .groupby("dmr_status")\
             .agg([
                 pl.count("beta").alias("count"),
                 pl.min("beta").alias("min"),
                 pl.max("beta").alias("max"),
                 pl.col("beta").quantile(0.25).alias("q25"),
                 pl.median("beta").alias("median"),
                 pl.col("beta").quantile(0.75).alias("q75"),
                 pl.mean("beta").alias("mean"),
                 pl.std("beta").alias("std"),
                 pl.var("beta").alias("variance"),
             ])\
             .collect()\
             .with_columns([pl.lit(tissue).alias("region"),
                            pl.lit(ancestry).alias("ancestry")])


def main():
    # Calculate average stats
    df = pl.DataFrame()
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        for LOCAL in [False, True]:
            df = df.vstack(load_dmrs(tissue, LOCAL))
    # Output results
    df.write_csv("ancestry_dmr_summary.tsv", sep="\t")
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
