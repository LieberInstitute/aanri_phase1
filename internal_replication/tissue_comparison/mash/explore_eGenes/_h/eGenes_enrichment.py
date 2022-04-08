"""
This script performs enrichment analysis on eGenes from mash model.
"""
import pandas as pd
from os import environ
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '32'
config = {
    'Caudate': "/ceph/users/jbenja13/github_projects/LieberInstituteBrainSeq" +\
    "Phase3CaudateSchizophrenia/analysis/eqtl/ancestry_comparison/" +\
    "meta_analysis/_m/genes/lfsr_allpairs_ancestry.txt.gz"
}

@lru_cache()
def get_eqtls(tissue):
    return pd.read_csv(config[tissue], sep='\t')


@lru_cache()
def get_eGenes(tissue):
    eqtls = get_eqtls(tissue)[(get_eqtls(tissue)["AA"] < 0.05) |
                              (get_eqtls(tissue)["EA"] < 0.05)].copy()
    return eqtls.sort_values(["AA", "EA"])\
                .groupby("gene_id").first().reset_index()


@lru_cache()
def get_mash_model():
    return pd.read_csv("../../_m/genes/lfsr_feature_4tissues.txt.gz", sep='\t')


@lru_cache()
def get_popDE(tissue):
    return get_mash_model()[(get_mash_model()[tissue] < 0.05)].copy()


def cal_enrichment(tissue):
    """
    Calculates Fisher's Exact test in universe u.
    Input tissue to subset data.
    """
    u = set(get_mash_model().Effect)
    a = set(get_eGenes(tissue).gene_id)
    b = set(get_popDE(tissue).Effect)
    yes_a = u.intersection(a); yes_b = u.intersection(b)
    no_a = u - a; no_b = u - b
    m = [[len(yes_a.intersection(yes_b)), len(no_a.intersection(yes_b))],
         [len(yes_a.intersection(no_b)), len(no_a.intersection(no_b))]]
    print("{:.1%} of popDE also eGenes in the {}!"\
          .format(len(yes_a.intersection(yes_b)) / len(yes_b), tissue))
    return fisher_exact(m)


def main():
    or_lt = []; pval_lt = []; tissue_lt = []
    for tissue in ["Caudate"]:
        oddratio, pvals = cal_enrichment(tissue)
        or_lt.append(oddratio); pval_lt.append(pvals); tissue_lt.append(tissue)
    fdr = multipletests(pval_lt, method='fdr_bh')[1]
    pd.DataFrame({"Tissue": tissue_lt, "OR": or_lt, "P-value": pval_lt,
                  "FDR": fdr})\
      .to_csv("popDE_eGene_enrichment.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
