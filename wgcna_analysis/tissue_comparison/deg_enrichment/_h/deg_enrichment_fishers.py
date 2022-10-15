# Examine enrichment in psychiatric disorders TWAS and DEGs
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

@lru_cache()
def get_background(tissue):
    fn = here("differential_analysis/tissue_comparison/summary_table/_m/",
              "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Tissue"] == tissue) & (df["Type"] == "Gene")].copy()


@lru_cache()
def get_ancestry_deg(tissue, direction):
    fn = here("differential_analysis/tissue_comparison/summary_table/_m/",
              "BrainSeq_ancestry_4features_4regions.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    df = df[(df["Tissue"] == tissue) & (df["Type"] == "Gene")].copy()
    df["ensemblID"] = df.gencodeID.str.replace("\\..*", "", regex=True)
    if direction=="all":
        return df
    elif direction=="up":
        return df[(df["posterior_mean"] > 0)].copy()
    else:
        return df[(df["posterior_mean"] < 0)].copy()


@lru_cache()
def get_wgcna_modules(tissue):
    if tissue == "Dentate Gyrus":
        new_tissue = "dentateGyrus"
    else:
        new_tissue = tissue.lower()
    fn = here("wgcna_analysis/aa_only/%s" % new_tissue,
              "_m/modules.csv")
    return pd.read_csv(fn, index_col=0)


def cal_enrichment(tissue, direction, mod):
    """
    Calculates Fisher's Exact test.
    Inputs: brainr region and DE direction of effect.
    """
    universe = set(get_background(tissue).Effect)
    de_set = set(get_ancestry_deg(tissue, direction).gencodeID)
    wgcna_df = get_wgcna_modules(tissue)
    wgcna_set = set(wgcna_df[(wgcna_df["module"] == mod)].index)
    yes_de = universe.intersection(de_set)
    yes_wgcna = universe.intersection(wgcna_set)
    no_de = universe - de_set
    no_wgcna = universe - wgcna_set
    m = [[len(yes_de.intersection(yes_wgcna)),
          len(no_de.intersection(yes_wgcna))],
         [len(yes_de.intersection(no_wgcna)),
          len(no_de.intersection(no_wgcna))]]
    return fisher_exact(m)


def enrichment_loop():
    dir_dict = {"all": "All", "up": "Downregulated in AA",
                "down": "Upregulated in AA"}
    dt = pd.DataFrame()
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        wgcna_df = get_wgcna_modules(tissue)
        modules = wgcna_df.module.unique()
        or_lt = []; pval_lt = []; dir_lt = []; mod_lt = []; mod_n = [];
        for mod in modules:
            wgcna_set = set(wgcna_df[(wgcna_df["module"] == mod)].index)
            for direction in ["all", "up", "down"]:
                oddratio, pvals = cal_enrichment(tissue, direction, mod)
                or_lt.append(oddratio); pval_lt.append(pvals);
                dir_lt.append(dir_dict[direction]); mod_lt.append(mod)
                mod_n.append(len(wgcna_set));
        fdr = multipletests(pval_lt, method='fdr_bh')[1]
        dtx = pd.DataFrame({"Tissue": tissue, "Module_ID": mod_lt,
                            "N_Genes":mod_n, "OR": or_lt, "P-value": pval_lt,
                            "FDR": fdr, "Direction":dir_lt})
        dt = pd.concat([dt, dtx], axis=0)
    return dt


def main():
    ## Direction dictionary
    dt = enrichment_loop()
    ## Save enrichment
    dt.to_csv("module_enrichment_analysis_DEGs.tsv",
              sep='\t', index=False)
    ## Session information
    session_info.show()


if __name__ == '__main__':
    main()
