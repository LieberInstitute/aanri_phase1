# Examine enrichment in psychiatric disorders TWAS and DEGs

import polars as pl
import session_info
from os import environ
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '10'

@lru_cache()
def get_ancestry_degs():
    fn = "../../../summary_table/_m/"+\
        "BrainSeq_ancestry_4features_4regions.txt.gz"
    return pl.read_csv(fn, separator="\t")\
             .filter(pl.col("Type") == "Gene")\
             .with_columns(pl.col("gencodeID")\
                           .str.replace("\\..*", "")\
                           .alias("ensembl_gene_id"))\
             .select(pl.col(["Tissue", "ensembl_gene_id", "Symbol", "lfsr"]),
                     pl.col("posterior_mean").alias("beta"))\
             .drop_nulls()


def load_data():
    fn = here("input/public/sun_microglia",
              "ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt")
    return pl.read_csv(fn, separator="\t")\
             .select([pl.col("cluster").apply(lambda x: f"MG{x}"),
                      pl.col("gene")])


def merge_data():
    return get_ancestry_degs().join(load_data(),
                                    left_on="Symbol",
                                    right_on="gene")


def main():
    ## Generate enrichment
    dt = merge_data()
    ## Save enrichment
    dt.write_csv("BrainSEQ.ancestryDEGs.activated_microglia.overlap.tsv",
                 separator='\t')
    ## Session information
    session_info.show()


if __name__ == '__main__':
    main()
