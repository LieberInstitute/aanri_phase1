"""
This script summarized DE results using mash model.
"""
import session_info
import polars as pl
from pyhere import here


def get_partial(tissue, feature):
    fn = f"../../../{tissue.lower()}/_m/{feature.lower()}s/partial_r2.tsv"
    return pl.scan_csv(fn, sep="\t", dtypes={"partial_r2": float})\
             .select([
                 pl.col("feature_id"), pl.col("partial_r2"),
                 pl.when(pl.col("partial_r2") < 0).then(0)\
                 .otherwise(pl.col("partial_r2")).alias("delta_partial")
             ])\
             .with_columns([pl.lit(tissue).alias("region"),
                            pl.lit(feature).alias("feature_type")])\
             .collect()


def merge_data():
    df = pl.DataFrame()
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        for feature in ["Gene", "Transcript"]:
            df = df.vstack(get_partial(tissue, feature))
    return df


def summarize_partial_r2():
    return merge_data().groupby(["region", "feature_type"])\
                       .agg([
                           pl.col("delta_partial").quantile(0.25).alias("q25"),
                           pl.median("delta_partial").alias("median"),
                           pl.col("delta_partial").quantile(0.75).alias("q75"),
                           pl.mean("delta_partial").alias("mean"),
                           pl.std("delta_partial").alias("std"),
                       ])\
                       .sort(["feature_type", "region"])


def get_mash():
    fn = here("differential_analysis/tissue_comparison/summary_table",
              "_m/BrainSeq_ancestry_4features_4regions.txt.gz")
    return pl.read_csv(fn, sep="\t")\
             .filter(pl.col("Type").is_in(["Gene", "Transcript"]))\
             .select([
                 pl.col("Tissue"), pl.col("Effect"),
                 pl.col("posterior_mean"), pl.col("Type"),
                 pl.when(pl.col("posterior_mean") > 0).then("EA bias").otherwise("AA bias").alias("direction")
             ])


def get_snp_assoc():
    fn = here("genetic_variation/elastic_net_resid5",
              "tissue_comparison/metrics_summary/_m",
              "genes.prediction_metrics_summary.txt.gz")
    return pl.read_csv(fn, sep="\t")


def extract_data(tissue):
    deg_df = get_mash().filter((pl.col("Type") == "Gene") &
                               (pl.col("Tissue") == tissue))\
                       .select([pl.all().exclude("Type")])
    vmr_df = merge_data().filter((pl.col("feature_type") == "Gene") &
                                 (pl.col("region") == tissue))\
                         .select([pl.all().exclude("feature_type")])
    snp_df = get_snp_assoc().filter(pl.col("region") == tissue)\
                            .select([pl.all().exclude("type")])
    return deg_df, snp_df, vmr_df


def feature_comparison(tissue):
    deg, snp, vmr = extract_data(tissue)
    snp_deg = set(deg.to_pandas().Effect) & set(snp.to_pandas().gene)
    snp_vmr = set(vmr.to_pandas().feature_id) & set(snp.to_pandas().gene)
    vmr_deg = set(vmr.to_pandas().feature_id) & set(deg.to_pandas().Effect)
    return None


def main():
    ## Save data for plotting
    merge_data().write_csv("partial_r2_de_summary.tsv", sep="\t")
    ## Summarize results
    print(summarize_partial_r2())
    print("Partial R2 >= 0.1")
    r2 = merge_data().filter(pl.col("delta_partial") >= 0.1)\
                     .groupby(["region", "feature_type"])\
                     .agg([pl.count("delta_partial")])\
                     .sort(["region", "feature_type"])
    print(r2)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
