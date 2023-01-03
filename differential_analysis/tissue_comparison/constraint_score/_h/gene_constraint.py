"""
Examines the gene constraint score and tests for 'enrichment' in
ancestry-related DEGs.
"""
import numpy as np
import pandas as pd
from pyhere import here
from functools import lru_cache
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import fisher_exact, ttest_ind, pearsonr

@lru_cache()
def load_constraint():
    fn = here("input/database/_m/gnomad/",
              "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
    return pd.read_csv(fn, sep="\t", compression="gzip")\
             .loc[:, ["gene", "pLI", "oe_lof_upper", "oe_lof_upper_bin", "p"]]\
             .dropna()


@lru_cache()
def get_annotation(feature="genes"):
    config = {
        "genes": here("input/text_files_counts/_m/caudate/gene_annotation.tsv"),
        "transcripts": here("input/text_files_counts/_m/caudate/tx_annotation.tsv"),
        "exons": here("input/text_files_counts/_m/caudate/exon_annotation.tsv"),
        "junctions": here("input/text_files_counts/_m/caudate/jxn_annotation.tsv"),
    }
    return pd.read_csv(config[feature], sep='\t')


@lru_cache()
def get_mash_degs():
    return pd.read_csv("../../_m/genes/lfsr_feature_4tissues.txt.gz", sep='\t')


@lru_cache()
def subset_tissue(tissue):
    return get_mash_degs().loc[:, ["Effect", tissue]]\
                          .rename(columns={tissue: "lfsr"})


@lru_cache()
def annotate_degs(tissue):
    return subset_tissue(tissue)\
        .merge(get_annotation(), left_on="Effect", right_on="names")\
        .drop(["names", "seqnames","start", "end", "Length"], axis=1)


@lru_cache()
def annot_effect_size(tissue):
    return pd.read_csv("../../_m/genes/posterior_mean_feature_4tissues.txt.gz",
                       sep='\t')\
             .loc[:, ["Effect", tissue]]\
             .rename(columns={tissue: "beta"})\
             .merge(annotate_degs(tissue), on="Effect").dropna()


@lru_cache()
def merge_data(tissue):
    err = 0.000001
    df = annot_effect_size(tissue)\
        .merge(load_constraint(), left_on="Symbol", right_on="gene")
    df["-log10(lfsr)"] = -np.log10(df.loc[:, "lfsr"] + err)
    df["Tissue"] = tissue.replace("Dentate.Gyrus", "Dentate Gyrus")
    return df


def old_cal_fishers(dfx, binx):
    """
    This function preforms enrichment examining constrained (<5) versus
    non-constrained (>5) using the LOEUF upper bound bins.

    Depletion  (OR < 1): DEGs are within less constrained genes.
    Enrichment (OR > 1): DEGs are within more constrained genes.
    """
    ## Low values indicate more constrained
    table = [[sum((dfx["lfsr"]<0.05) & ((dfx["oe_lof_upper_bin"]==binx))),
              sum((dfx["lfsr"]<0.05) & ((dfx["oe_lof_upper_bin"]!=binx)))],
             [sum((dfx["lfsr"]>=0.05) & ((dfx["oe_lof_upper_bin"]==binx))),
              sum((dfx["lfsr"]>=0.05) & ((dfx["oe_lof_upper_bin"]!=binx)))]]
    #print(table)
    return fisher_exact(table)


def cal_fishers(dfx):
    """
    This function preforms enrichment examining constrained (<5) versus
    non-constrained (>5) using the LOEUF upper bound bins.

    Depletion  (OR < 1): DEGs are within less constrained genes.
    Enrichment (OR > 1): DEGs are within more constrained genes.
    """
    ## High values (>.9) indicate more constrained
    table = [[sum((dfx["lfsr"]<0.05) & ((dfx["pLI"]<0.9))),
              sum((dfx["lfsr"]<0.05) & ((dfx["pLI"]>=0.9)))],
             [sum((dfx["lfsr"]>=0.05) & ((dfx["pLI"]<0.9))),
              sum((dfx["lfsr"]>=0.05) & ((dfx["pLI"]>=0.9)))]]
    #print(table)
    return fisher_exact(table)


def constraint_comparison():
    ## This function does not take direction of effect to account.
    with open("statistical.log", "w") as f:
        print("Preform statistical analysis!", file=f)
        for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
            print("%s" % tissue, file=f)
            df = merge_data(tissue)
            tstat, pval1 = ttest_ind(df.loc[:, "lfsr"], df.loc[:, "oe_lof_upper"])
            print("T-test:\t\t T-stat > %.3f, p-value < %.1e" %(tstat, pval1), file=f)
            rho, pval2 = pearsonr(df.loc[:, "lfsr"], df.loc[:, "oe_lof_upper"])
            print("Pearson corr: rho > %.3f, p-value < %.1e" % (rho, pval2), file=f)        


def old_enrichment_analysis():
    ## Enrichment analysis
    dt = pd.DataFrame(); dft = pd.DataFrame()
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        df = merge_data(tissue)
        tissue_lt = []; oddratio = []; pvalues = []; bins = []
        for upper_bin in range(0,10):
            odds, pval = old_cal_fishers(df, upper_bin)
            bins.append(upper_bin); oddratio.append(odds)
            pvalues.append(pval); tissue_lt.append(tissue)
        fdr = fdrcorrection(pvalues)[1]
        dt = pd.concat([dt, pd.DataFrame({"Tissue": tissue_lt, "Upper_Bin": bins,
                                          "Odds_Ratio": oddratio, "P_value": pvalues,
                                          "FDR": fdr})])
        dft = pd.concat([dft, df])
    return dt, dft


def enrichment_analysis():
    ## Enrichment analysis
    dft = pd.DataFrame(); tissue_lt = []; oddratio = []; pvalues = []
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        df = merge_data(tissue)
        odds, pval = cal_fishers(df)
        oddratio.append(odds)
        pvalues.append(pval)
        tissue_lt.append(tissue)
        dft = pd.concat([dft, df])
    fdr = fdrcorrection(pvalues)[1]
    dt = pd.DataFrame({"Tissue": tissue_lt, "Odds_Ratio": oddratio,
                       "P_value": pvalues, "FDR": fdr})
    return dt, dft


def main():
    dt, dft = enrichment_analysis()
    dx, _ = old_enrichment_analysis()
    dt.to_csv("constrain_enrichment_pLI.tsv", sep='\t', index=False)
    dx.to_csv("constrain_enrichment.tsv", sep='\t', index=False)
    dft.to_csv("brainseq_degs_constrain_score.tsv", sep='\t', index=False)
    ## Additional statistics
    constraint_comparison()
    

if __name__ == '__main__':
    main()
