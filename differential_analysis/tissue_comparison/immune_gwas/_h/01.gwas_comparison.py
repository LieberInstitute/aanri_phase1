#### This script cleans immune GWAS data for
#### public comparison.

import session_info
import polars as pl
from pyhere import here
from functools import lru_cache

@lru_cache()
def get_ancestry_degs():
    fn = "../../summary_table/_m/"+\
        "BrainSeq_ancestry_4features_4regions.txt.gz"
    return pl.read_csv(fn, separator="\t")\
             .filter(pl.col("Type") == "Gene")


def subset_region(region):
    return get_ancestry_degs()\
        .filter(pl.col("Tissue") == region)\
        .select(pl.col(["gencodeID", "Symbol", "lfsr"]),
                pl.col("posterior_mean").alias("beta"))\
        .drop_nulls()


def get_mhc():
    fn = "../../mhc_enrichment/_m/mhc_genes.csv"
    return pl.read_csv(fn)\
        .select(pl.col("gene_id", "gene_name"))\
        .with_columns(pl.lit(1).alias("mhc"))


def load_data():
    fn1 = here("input/immune_gwas/_m", "orru_2013_gwas.csv")
    df1 = pl.read_csv(fn1)\
            .select(
                pl.col("Trait"),
                pl.col(" Candidate genes")\
                .str.split(by=",")\
                .alias("candidate_genes")\
                .list.to_struct()
            ).unnest("candidate_genes")
    fn2 = here("input/immune_gwas/_m", "patin_2018_gwas.csv")
    df2 = pl.scan_csv(fn2)\
            .select(pl.col("Immunophenotype"),
                    pl.col("Candidate gene").alias("candidate_genes"))\
            .collect()
    fn3 = here("input/immune_gwas/_m", "orru_2020_gwas.tsv")
    df3 = pl.scan_csv(fn3, separator="\t",
                      infer_schema_length=10000)\
            .select(pl.col(["DISEASE_TRAIT", "INITIAL_SAMPLE_SIZE",
                            "CHR_ID", "CHR_POS", "MAPPED_GENE",
                            "SNP_GENE_IDS"]))\
            .collect().drop_nulls()\
            .select(pl.col("DISEASE_TRAIT"),
                    pl.col("MAPPED_GENE")\
                    .str.split(by=", ")\
                    .alias("candidate_genes")\
                    .list.to_struct())\
            .unnest("candidate_genes")
    return df1, df2, df3


def gwas_immune():
    ## Load GWAS
    df1, df2, df3 = load_data()
    ## Clean GWAS
    genes1 = df1.select(pl.col("Trait").alias("trait"),
                        pl.col("field_0").alias("candidate_genes"))\
                .with_columns(pl.lit("Orru_2013").alias("dataset"))
    genes2 = df2.select(pl.col("Immunophenotype").alias("trait"),
                        pl.col("candidate_genes"))\
                .with_columns(pl.lit("Patin_2018").alias("dataset"))
    genes3 = df3.select(pl.col("DISEASE_TRAIT").alias("trait"),
                        pl.col("field_0").alias("candidate_genes"))\
                .with_columns(pl.lit("Orru_2020").alias("dataset"))
    return pl.concat([genes1, genes2, genes3], how="vertical")


def comparison_immune(region, MHC):
    ## Load data
    df1, df2, df3 = load_data()
    dx = subset_region(region)\
        .join(get_mhc(), right_on="gene_id",
              left_on="gencodeID", how="left")\
        .with_columns(pl.when(pl.col("mhc") == None)\
                      .then(0).otherwise(pl.col("mhc"))\
                      .keep_name())
    if MHC:
        ## Select only MHC region
        dx = dx.filter(pl.col("mhc") == 1).to_pandas()
    else:
        ## Exclude MHC region
        dx = dx.filter(pl.col("mhc") == 0).to_pandas()
    ## Overlap between immune GWAS    
    genes1 = df1.filter(pl.col("field_0").is_in(dx.Symbol))\
                .select(pl.col("Trait").alias("trait"),
                        pl.col("field_0").alias("candidate_genes"))\
                .with_columns(pl.when(MHC).then("Yes").otherwise("No").alias("mhc"))
    genes2 = df2.filter(pl.col("candidate_genes").is_in(dx.Symbol))\
                .select(pl.col("Immunophenotype").alias("trait"),
                        pl.col("candidate_genes"))\
                .with_columns(pl.when(MHC).then("Yes").otherwise("No").alias("mhc"))
    genes3 = df3.filter(pl.col("field_0").is_in(dx.Symbol))\
                .select(pl.col("DISEASE_TRAIT").alias("trait"),
                        pl.col("field_0").alias("candidate_genes"))\
                .with_columns(pl.when(MHC).then("Yes").otherwise("No").alias("mhc"))
    return pl.concat([genes1, genes2, genes3], how="vertical")\
        .with_columns(pl.lit(region).alias("region"))


def main():
    ## Combine analysis
    df = pl.DataFrame()
    for region in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        for MHC in [True, False]:
            df = pl.concat([df, comparison_immune(region, MHC)],
                           how="vertical")
    ## Save results
    df.write_csv("immune_gwas.overlap.ancestry_degs.tsv", separator="\t")
    gwas_immune().write_csv("immune_gwas.annotation.tsv", separator="\t")
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
