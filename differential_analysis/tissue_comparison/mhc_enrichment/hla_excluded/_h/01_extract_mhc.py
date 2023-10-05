#### This script extracts the MHC region based on hg38 annotation

import session_info
import polars as pl
import pyranges as pr
from pyhere import here
from gtfparse import read_gtf


def extract_region():
    return pr.PyRanges(chromosomes = ['chr6'],
                       starts=[25_726_063],
                       ends=[33_400_644])


def get_gtf():
    fn = here("input/genomes/_m/gencode.v25.annotation.gtf")
    df = read_gtf(str(fn))
    return df.filter(pl.col("feature") == "gene")\
             .select([
                 pl.col("seqname").alias("Chromosome"),
                 pl.col("start").alias("Start"),
                 pl.col("end").alias("End"),
                 pl.col("gene_id"), pl.col("gene_name"),
             ])


def main():
    # Main
    pr_genes = pr.PyRanges(get_gtf().to_pandas())
    pr_mhc = extract_region()
    xmhc_genes = pr_genes.overlap(pr_mhc, strandedness=False,how='first')\
                         .as_df().rename(columns=str.lower)
    hla_genes = xmhc_genes[xmhc_genes.gene_name.str.contains("HLA")].copy()
    hla_genes.to_csv("hla_genes.csv", index=False)
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()

    
