# Examine enrichment in psychiatric disorders TWAS and DEGs
import pandas as pd
import session_info
from os import environ
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '10'

@lru_cache()
def get_eqtl_interacting(tissue):
    new_tissue = tissue.replace(" ", "_")
    fn = "../../_m/lfsr_allpairs_ancestry.txt.gz"
    df0 = pd.read_csv(fn, sep='\t', nrows=100,
                      usecols=["gene_id", new_tissue])
    return pd.read_csv(fn, sep='\t', dtype=df0.dtypes.to_dict(),
                       usecols=["gene_id", new_tissue],
                       compression="gzip")


@lru_cache()
def get_eqtl_main(tissue):
    new_tissue = tissue.replace(" ", "_")
    fn = "../../../../tissue_comparison/transcripts/"+\
        "_m/lfsr_allpairs_ancestry.txt.gz"
    df0 = pd.read_csv(fn, sep='\t', nrows=100,
                      usecols=["gene_id", new_tissue])
    return pd.read_csv(fn, sep='\t', dtype=df0.dtypes.to_dict(),
                       usecols=["gene_id", new_tissue],
                       compression="gzip")


@lru_cache()
def get_background(tissue):
    eqtl_df1 = get_eqtl_interacting(tissue)
    eqtl_df2 = get_eqtl_main(tissue)
    return list(set(eqtl_df1.gene_id) | set(eqtl_df2.gene_id))


@lru_cache()
def get_eGene_interacting(tissue):
    new_tissue = tissue.replace(" ", "_")
    return get_eqtl_interacting(tissue)[(get_eqtl_interacting(tissue)[new_tissue] < 0.05)]\
        .drop_duplicates(subset="gene_id")


@lru_cache()
def get_eGene_main(tissue):
    new_tissue = tissue.replace(" ", "_")
    return get_eqtl_main(tissue)[(get_eqtl_main(tissue)[new_tissue] < 0.05)]\
        .drop_duplicates(subset="gene_id")


def cal_enrichment(tissue):
    """
    Calculates Fisher's Exact test.
    Inputs: brainr region and DE direction of effect.
    """
    universe = set(get_background(tissue))
    eGene_main = set(get_eGene_main(tissue).gene_id);
    eGene_interact = set(get_eGene_interacting(tissue).gene_id)
    yes_main = universe.intersection(eGene_main)
    yes_interact = universe.intersection(eGene_interact)
    no_main = universe - eGene_main
    no_interact = universe - eGene_interact
    m = [[len(yes_main.intersection(yes_interact)),
          len(no_main.intersection(yes_interact))],
         [len(yes_main.intersection(no_interact)),
          len(no_main.intersection(no_interact))]]
    return fisher_exact(m)


def enrichment_loop():
    dt = pd.DataFrame()
    or_lt = []; pval_lt = []; tissue_lt = []
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        oddratio, pvals = cal_enrichment(tissue)
        or_lt.append(oddratio); pval_lt.append(pvals);
        tissue_lt.append(tissue)
    fdr = multipletests(pval_lt, method='fdr_bh')[1]
    return pd.DataFrame({"Tissue": tissue_lt, "OR": or_lt,
                         "P-value": pval_lt, "FDR": fdr})


def main():
    ## Enrichment analysis
    dt = enrichment_loop()
    ## Save enrichment
    dt.to_csv("eGene_enrichment_fishers.tsv", sep='\t', index=False)
    session_info.show()


if __name__ == '__main__':
    main()
