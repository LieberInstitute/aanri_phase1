# Examine absolute allele frequency differences between DEGs and other eGenes
import pandas as pd
import session_info
from os import environ
from pyhere import here
from functools import lru_cache
from scipy.stats import mannwhitneyu

environ['NUMEXPR_MAX_THREADS'] = '10'

@lru_cache()
def get_eqtl(tissue):
    new_tissue = tissue.replace(" ", "_")
    fn = here("eqtl_analysis/tissue_comparison/genes",
              "_m/lfsr_allpairs_ancestry.txt.gz")
    df0 = pd.read_csv(fn, sep='\t', nrows=100,
                      usecols=["gene_id", "variant_id", new_tissue])
    return pd.read_csv(fn, sep='\t', dtype=df0.dtypes.to_dict(),
                       usecols=["gene_id", "variant_id", new_tissue],
                       compression="gzip")


@lru_cache()
def get_eGene(tissue, fdr=1):
    new_tissue = tissue.replace(" ", "_")
    eqtl_df = get_eqtl(tissue)[(get_eqtl(tissue)[new_tissue] < fdr)].copy()
    snp_ids = eqtl_df.variant_id.str.split("_", expand=True)
    eqtl_df["rsid"] = snp_ids[4]
    return eqtl_df.rename(columns={new_tissue:"lfsr"})


@lru_cache()
def top_eqtl(tissue):
    return get_eGene(tissue).sort_values("lfsr")\
                            .groupby("gene_id")\
                            .first().reset_index()


@lru_cache()
def get_ancestry_deg(tissue, fdr=1):
    fn = here("differential_analysis/tissue_comparison/summary_table/",
              "_m/BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    df = df[(df["Tissue"] == tissue) & (df["Type"] == "Gene") &
            (df["lfsr"] < fdr)].copy()
    return df.loc[:, ["Effect", "lfsr", "posterior_mean"]].copy()


@lru_cache()
def get_allele_freq():
    df = pd.DataFrame()
    for chrom in range(1, 23):
        fn = f"../_m/allele_freq_difference.1000GP.chr{chrom}.tsv"
        df = pd.concat([df, pd.read_csv(fn, sep='\t', usecols=[1,4])], axis=0)
    return df

    
def allele_freq_topSNP(tissue):
    new_tissue = tissue.lower().replace(" ", "_")
    eqtl_df = top_eqtl(tissue).merge(get_allele_freq(),
                                     left_on="rsid", right_on="SNP")\
                              .merge(get_ancestry_deg(tissue),
                                     left_on="gene_id", right_on="Effect",
                                     suffixes=('_eQTL', '_DE'))
    degs = eqtl_df[(eqtl_df["lfsr_DE"] < 0.05)].copy()
    others = eqtl_df[(eqtl_df["lfsr_DE"] > 0.05)].copy()
    df1 = degs; df1["DEG"] = 1; df2 = others; df2["DEG"] = 0
    pd.concat([df1, df2], axis=0)\
      .to_csv(f"{new_tissue}.eSNP.degs_AFD.tsv", sep="\t", index=False)
    stat, pval = mannwhitneyu(degs.AFD, others.AFD, alternative="greater")
    mu1 = degs.AFD.mean(); mu2 = others.AFD.mean()
    sd1 = degs.AFD.std(); sd2 = others.AFD.std()
    return pd.DataFrame({"Feature": ["eSNP"], "Pvalue": [pval],
                         "DEG_mean": [mu1], "DEG_sd": [sd1],
                         "Other_mean": [mu2], "Other_sd": [sd2]})


def allele_freq_eGene(tissue):
    new_tissue = tissue.lower().replace(" ", "_")
    eqtl_df = get_eGene(tissue, 0.05).merge(get_allele_freq(),
                                            left_on="rsid", right_on="SNP")
    degs = eqtl_df[(eqtl_df["gene_id"].isin(get_ancestry_deg(tissue, 0.05).Effect))]\
        .loc[:, ["gene_id", "AFD"]].groupby("gene_id").mean().reset_index()
    # degs["AFD"] = degs.AFD / 2
    others = eqtl_df[~(eqtl_df["gene_id"].isin(get_ancestry_deg(tissue, 0.05).Effect))]\
        .loc[:, ["gene_id", "AFD"]].groupby("gene_id").mean().reset_index()
    df1 = degs; df1["DEG"] = 1; df2 = others; df2["DEG"] = 0
    pd.concat([df1, df2], axis=0)\
      .to_csv(f"{new_tissue}.eGene.degs_AFD.tsv", sep="\t", index=False)
    # others["AFD"] = others.AFD / 2
    stat, pval = mannwhitneyu(degs.AFD, others.AFD, alternative="greater")
    mu1 = degs.AFD.mean(); mu2 = others.AFD.mean()
    sd1 = degs.AFD.std(); sd2 = others.AFD.std()
    return pd.DataFrame({"Feature": ["eGene"], "Pvalue": [pval],
                         "DEG_mean": [mu1], "DEG_sd": [sd1],
                         "Other_mean": [mu2], "Other_sd": [sd2]})


def main():
    ## Calculate differences
    dt = pd.DataFrame()
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        df = pd.concat([allele_freq_eGene(tissue),
                        allele_freq_topSNP(tissue)], axis=0)
        df["Region"] = tissue
        dt = pd.concat([dt, df], axis=0)
    ## Save enrichment
    dt.to_csv("allele_freq.DEGs_Others.tsv", sep='\t', index=False)
    session_info.show()


if __name__ == '__main__':
    main()
