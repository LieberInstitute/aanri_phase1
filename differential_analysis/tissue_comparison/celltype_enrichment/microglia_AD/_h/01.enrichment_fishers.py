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
    fn = here("input/public/sun_microglia",
              "ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt")
    return pl.read_csv(fn, separator="\t")\
             .select([pl.col("cluster"), pl.col("gene")])


def cal_enrichment(tissue, direction, cluster):
    """
    Calculates Fisher's Exact test.
    """
    universe = set(get_database().external_gene_name)
    ancestry = set(subset_degs(tissue, direction).get_column("Symbol"))
    ct = set(load_data().filter(pl.col("cluster") == cluster)\
                 .get_column("gene"))
    yes_a = universe.intersection(ancestry); yes_i = universe.intersection(ct)
    no_a = universe - ancestry; no_i = universe - ct
    m = [[len(yes_a.intersection(yes_i)), len(no_a.intersection(yes_i))],
         [len(yes_a.intersection(no_i)), len(no_a.intersection(no_i))]]
    return fisher_exact(m)


def enrichment_loop():
    dir_dict = {"all": "All DEGs", "up": "Increased in EA",
                "down": "Increased in AA"}
    clusters = load_data().get_column("cluster").unique()
    dt = pl.DataFrame()
    for direction in ["all", "up", "down"]:
        or_lt = []; pval_lt = []; tissue_lt = []; data_lt = []
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            for cluster in clusters:
                oddratio, pvals = cal_enrichment(tissue, direction, cluster)
                or_lt.append(oddratio); pval_lt.append(pvals);
                tissue_lt.append(tissue); data_lt.append(cluster)
        fdr = multipletests(pval_lt, method='fdr_bh')[1]
        dtx = pl.DataFrame({"Tissue": tissue_lt, "Cluster": data_lt,
                            "OR": or_lt, "P-value": pval_lt, "FDR": fdr,
                            "Direction":dir_dict[direction]})
        dt = pl.concat([dt, dtx], how="vertical")
    return dt


def main():
    ## Generate enrichment
    dt = enrichment_loop()
    ## Save enrichment
    dt.write_csv("celltype_enrichment.tsv", separator='\t')
    ## Session information
    session_info.show()


if __name__ == '__main__':
    main()
