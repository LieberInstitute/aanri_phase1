"""
This script preforms the Fisher's pairwise enrichment analysis
between brain regions.
"""
import pandas as pd
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

@lru_cache()
def get_public_deg():
    fn = "../../../../input/public/nedelec_immune_bulk.xlsx"
    return pd.read_excel(fn, sheet_name="B) pop-DE", skiprows=3)


@lru_cache()
def get_mash():
    # get effect size
    fn = "../../_m/genes/lfsr_feature_4tissues.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    df["ensembl_gene_id"] = df.Effect.str.replace("\\..*", "", regex=True)
    return df


@lru_cache()
def merge_data():
    """
    Merge data to get shared features to check for enrichment.
    """
    return pd.merge(get_mash(), get_public_deg(), on="ensembl_gene_id")


def cal_fishers(tissue, lab):
    table = [[sum((merge_data()[tissue]<0.05) & ((merge_data()[lab]<0.05))),
              sum((merge_data()[tissue]<0.05) & ((merge_data()[lab]>=0.05)))],
             [sum((merge_data()[tissue]>=0.05) & ((merge_data()[lab]<0.05))),
              sum((merge_data()[tissue]>=0.05) & ((merge_data()[lab]>=0.05)))]]
    print(table)
    return fisher_exact(table, alternative='greater')


def enrichment_loop():
    comp_list = {"NI_fdr": "Non-infected", "L_fdr": "Listeria",
                 "S_fdr": "Salmonella"}
    or_lt = []; pval_lt = []; tissue_lt = []; comp_lt = []
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        for comp in ["NI_fdr", "L_fdr", "S_fdr"]:
            ## Correlated effect sizes
            odds, pval = cal_fishers(tissue, comp)
            or_lt.append(odds); pval_lt.append(pval);
            tissue_lt.append(tissue); comp_lt.append(comp_list[comp])
    fdr = multipletests(pval_lt, method='fdr_bh')[1]
    dt  = pd.DataFrame({"Tissue": tissue_lt, "Conditions": comp_lt,
                        "OR": or_lt, "P-value": pval_lt, "FDR": fdr})
    return dt


def main():
    ## Enrichment loop and save results
    enrichment_loop()\
        .to_csv("nedelec_immune_enrichment_analysis_DEGs.tsv",
                sep='\t', index=False)

    ## Calculate rho
    with open("fishers_enrichment.log", "w") as f:
        for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
            for label in ["NI_fdr", "L_fdr", "S_fdr"]:
                ## Correlated effect sizes
                odds, pval = cal_fishers(tissue, label)
                print("%s vs %s:\t  Odds Ratio > %.3f, p-value < %.1e" %
                      (tissue, label, odds, pval), file=f)


if __name__ == "__main__":
    main()
