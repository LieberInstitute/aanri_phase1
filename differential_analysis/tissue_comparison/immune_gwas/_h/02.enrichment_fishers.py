# Examine enrichment in psychiatric disorders TWAS and DEGs

import polars as pl
import session_info
from os import environ
from pyhere import here
from pybiomart import Dataset
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '10'

@lru_cache()
def get_database():
    dataset = Dataset(name="hsapiens_gene_ensembl",
                      host="http://www.ensembl.org",
                      use_cache=True)
    db = dataset.query(attributes=["ensembl_gene_id",
                                   "external_gene_name",
                                   "entrezgene_id"],
                       use_attr_names=True).dropna(subset=['entrezgene_id'])
    return db


@lru_cache()
def get_ancestry_degs():
    fn = "../../summary_table/_m/"+\
        "BrainSeq_ancestry_4features_4regions.txt.gz"
    return pl.read_csv(fn, separator="\t")\
             .filter(pl.col("Type") == "Gene")\
             .with_columns(pl.col("gencodeID")\
                           .str.replace("\\..*", "")\
                           .alias("ensembl_gene_id"))\
             .select(pl.col(["Tissue", "ensembl_gene_id", "Symbol", "lfsr"]),
                     pl.col("posterior_mean").alias("beta"))\
             .drop_nulls()


@lru_cache()
def subset_degs(tissue, direction):
    if direction=="all":
        return get_ancestry_degs()\
            .filter(pl.col("Tissue") == tissue)
    elif direction=="up":
        return get_ancestry_degs()\
            .filter((pl.col("Tissue") == tissue) & (pl.col("beta") > 0))
    else:
        return get_ancestry_degs()\
            .filter((pl.col("Tissue") == tissue) & (pl.col("beta") < 0))


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


def cal_enrichment(tissue, direction, dataset):
    """
    Calculates Fisher's Exact test.
    """
    universe = set(get_database().external_gene_name)
    ancestry = set(subset_degs(tissue, direction).get_column("Symbol"))
    immune = set(gwas_immune().filter(pl.col("dataset") == dataset)\
                 .get_column("candidate_genes"))
    yes_a = universe.intersection(ancestry); yes_i = universe.intersection(immune)
    no_a = universe - ancestry; no_i = universe - immune
    m = [[len(yes_a.intersection(yes_i)), len(no_a.intersection(yes_i))],
         [len(yes_a.intersection(no_i)), len(no_a.intersection(no_i))]]
    return fisher_exact(m)


def enrichment_loop():
    dir_dict = {"all": "All", "up": "Increased in EA",
                "down": "Increased in AA"}
    datasets = gwas_immune().get_column("dataset").unique()
    dt = pl.DataFrame()
    for direction in ["all", "up", "down"]:
        or_lt = []; pval_lt = []; tissue_lt = []; data_lt = []
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            for dataset in datasets:
                oddratio, pvals = cal_enrichment(tissue, direction, dataset)
                or_lt.append(oddratio); pval_lt.append(pvals);
                tissue_lt.append(tissue); data_lt.append(dataset)
        fdr = multipletests(pval_lt, method='fdr_bh')[1]
        dtx = pl.DataFrame({"Tissue": tissue_lt, "Dataset": data_lt,
                            "OR": or_lt, "P-value": pval_lt, "FDR": fdr,
                            "Direction":dir_dict[direction]})
        dt = pl.concat([dt, dtx], how="vertical")
    return dt


def main():
    ## Generate enrichment
    dt = enrichment_loop()
    ## Save enrichment
    dt.write_csv("immune_GWAS_enrichment.tsv", separator='\t')
    ## Session information
    session_info.show()


if __name__ == '__main__':
    main()
