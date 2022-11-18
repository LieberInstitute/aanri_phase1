## Examine enrichment in psychiatric disorders TWAS and DE transcripts
## Using a fitted regularized maximum likelihood model.

import numpy as np
import pandas as pd
import session_info
from os import environ
import statsmodels.api as sm
from functools import lru_cache
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '10'
config = {
    "caud8_file": "/dcs04/lieber/statsgen/jbenjami/resources/degs/" +\
    "brainseq_phase3/BrainSeq_Phase3_Caudate_DifferentialExpression_DxSZ_all.txt.gz",
    "dlpfc_file": "/dcs04/lieber/statsgen/jbenjami/resources/degs/" +\
    "brainseq_phase2/Brainseq_phase2_qsvs_age17_noHGold_sz_DEG.tsv",
    "hippo_file": "/dcs04/lieber/statsgen/jbenjami/resources/degs/" +\
    "brainseq_phase2/Brainseq_phase2_qsvs_age17_noHGold_sz_DEG.tsv",
    'cmc_file': '/dcs04/lieber/statsgen/jbenjami/resources/degs/' +\
    'cmc/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500'+\
    '_gene-adjustedSVA-differentialExpression-includeAncestry-DxSCZ-DE.tsv',
    'PW_parkinson_file': '/dcs04/lieber/statsgen/jbenjami/resources/degs/' +\
    'Nido2020_PD_pwCohort.csv',
    'alzheimers_file': '/dcs04/lieber/statsgen/jbenjami/resources/degs/' +\
    'MarquesCoelho2021_AD.csv',
    "libd_bipolar_file": "/dcs04/lieber/statsgen/jbenjami/resources/degs/" +\
    "bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation_all.csv",
    'gandal_de_file': "/dcs04/lieber/statsgen/jbenjami/psychENCODE/expression_results/" +\
    "_m/gandal2018_psychENCODE_DE_results.xlsx",
    'twas_asd_file': '/dcs04/lieber/statsgen/jbenjami/psychENCODE/_m/psychENCODE_twas_asd.csv',
    'twas_sz_file': '/dcs04/lieber/statsgen/jbenjami/psychENCODE/_m/psychENCODE_twas_sz.csv',
    'twas_bd_file': '/dcs04/lieber/statsgen/jbenjami/psychENCODE/_m/psychENCODE_twas_bd.csv',
    'twas_bs_file': '/dcs04/lieber/statsgen/jbenjami/resources/twas/' +\
    'TWAS_gene_tissue_summary.csv',
    'twas_ad_file': '/dcs04/lieber/statsgen/jbenjami/resources/twas/Gockley2021_AD.tsv',
    'twas_bd_libd_file': '/dcs04/lieber/statsgen/jbenjami/resources/twas/' +\
    'twas_2regions_bpd_output_allTested.csv',
}

@lru_cache()
def get_deg(fn, REGION):
    df = pd.read_csv(fn, sep='\t')
    if REGION.lower() != "cmc":
        df = df[(df["Type"] == "Gene")].copy()
    if (REGION.lower() == "caudate") | (REGION.lower() == "cmc"):
        return df
    else:
        return df[(df["Tissue"] == REGION)].copy()


@lru_cache()
def get_sig_deg(fn, REGION):
    return get_deg(fn, REGION)[(get_deg(fn, REGION)["adj.P.Val"] < 0.05)].copy()


@lru_cache()
def get_bd_deg():
    df = pd.read_csv(config["libd_bipolar_file"], index_col=0)
    df["ensemblID"] = df.gencodeID.str.replace("\..*", "", regex=True)
    return df[(df["Type"] == "Gene")].copy()


@lru_cache()
def get_bd_sig_deg(label):
    return get_bd_deg()[(get_bd_deg()[label] < 0.05)].copy()


@lru_cache()
def get_gandal_deg():
    return pd.read_excel(config["gandal_de_file"], sheet_name="DGE")


@lru_cache()
def get_parkinson_deg():
    return pd.read_csv(config['PW_parkinson_file'])


@lru_cache()
def get_alzheimers_deg(dataset):
    df = pd.read_csv(config['alzheimers_file'], skiprows=1)\
           .drop_duplicates(subset="GENEID")
    df["ensemblID"] = df.GENEID.str.replace("\..*", "", regex=True)
    return df[(df["gene.padj"] < 0.05) & (df["dataset"] == dataset)].copy()


@lru_cache()
def get_ancestry_deg(tissue, direction):
    fn = "../../summary_table/_m/"+\
        "BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    df = df[(df["Tissue"] == tissue) &
            (df["Type"] == "Transcript")].copy()
    n = df.groupby("gencodeID").size().reset_index().rename(columns={0:"n"})
    df = df.merge(n, on="gencodeID")
    df["ensemblID"] = df.gencodeID.str.replace("\\..*", "", regex=True)
    df["popDE"] = [1 if x <= 0.05 else 0 for x in df.lfsr]
    if direction=="all":
        return df.loc[:, ["ensemblID","lfsr", "popDE", "n"]]
    elif direction=="up":
        return df[(df["posterior_mean"] > 0)]\
            .loc[:,["ensemblID","lfsr", "popDE", "n"]].copy()
    else:
        return df[(df["posterior_mean"] < 0)]\
            .loc[:,["ensemblID","lfsr", "popDE", "n"]].copy()


@lru_cache()
def get_bs_twas():
    return pd.read_csv(config["twas_bs_file"])


@lru_cache()
def get_libd_bd_twas():
    return pd.read_csv(config['twas_bd_libd_file'])


@lru_cache()
def get_gandal_twas(fn):
    df = pd.read_csv(fn)
    return df[(df["TWAS.Bonferroni"] < 0.05)].copy()


@lru_cache()
def get_ad_twas():
    df = pd.read_csv(config['twas_ad_file'], sep='\t', header=None,
                     names=["chr", "set", "feature", "start", "end", "na1",
                            "strand", "na2", "geneid", "gene_name"])
    df["ensemblID"] = df.geneid.str.replace("\..*", "", regex=True)
    return df


def extract_oddsratio(res):
    conf = res.conf_int()
    conf.columns = ['5%', '95%']
    conf['Odds Ratio'] = res.params
    return np.exp(conf)


def logreg_enrichment(dt, genes):
    """
    Calculates enrichment using logistic regression.
    """
    dx = dt.merge(pd.DataFrame({"ensemblID":genes, "CC":1}),
                  on="ensemblID",how="left").fillna(0)
    y = dx.CC
    x = dx.loc[:, ["popDE", "n"]]
    x = sm.add_constant(x)
    result = sm.Logit(y,x).fit_regularized()
    cil, ciu, oddratio = extract_oddsratio(result).loc["popDE"]
    pval = result.pvalues.loc["popDE"]
    return cil, ciu, oddratio, pval


def load_ancestry_de(direction="all"):
    caudate = get_ancestry_deg("Caudate", direction)
    gyrus   = get_ancestry_deg("Dentate Gyrus", direction)
    dlpfc   = get_ancestry_deg("DLPFC", direction)
    hippo   = get_ancestry_deg("Hippocampus", direction)
    return {"Caudate": caudate, "Dentate Gyrus": gyrus,
            "DLPFC": dlpfc, "Hippocampus": hippo}


def get_public_degs():
    # DEGs
    ## BrainSeq SZ
    bs_caudate_degs = list(get_sig_deg(config["caud8_file"], "Caudate").ensemblID)
    bs_dlpfc_degs = list(get_sig_deg(config["dlpfc_file"], "DLPFC").ensemblID)
    bs_hippo_degs = list(get_sig_deg(config["hippo_file"], "Hippocampus").ensemblID)
    ## CommonMind SZ
    cmc_dlpfc_degs = list(get_sig_deg(config["cmc_file"], "CMC").genes)
    ## LIBD Bipolar
    bs_sacc_degs = list(get_bd_sig_deg("adj.P.Val_sACC").ensemblID)
    bs_amyg_degs = list(get_bd_sig_deg("adj.P.Val_Amyg").ensemblID)
    ## PsychENCODE (Gandal)
    psy_sz = list(get_gandal_deg()[(get_gandal_deg()["SCZ.fdr"] < 0.05)]\
                  .ensembl_gene_id)
    psy_asd = list(get_gandal_deg()[(get_gandal_deg()["ASD.fdr"] < 0.05)]\
                   .ensembl_gene_id)
    psy_bd = list(get_gandal_deg()[(get_gandal_deg()["BD.fdr"] < 0.05)]\
                  .ensembl_gene_id)
    ## Parkinson's
    pd_degs = list(get_parkinson_deg()[(get_parkinson_deg()["padj"] < 0.05)]\
                   .EnsemblID)
    ## Alzheimer's
    ad_mayo = list(get_alzheimers_deg("MAYO").ensemblID)
    ad_rosmap = list(get_alzheimers_deg("ROSMAP").ensemblID)
    ad_msbb_mb36 = list(get_alzheimers_deg("MSBB BM36").ensemblID)
    ad_msbb_mb44 = list(get_alzheimers_deg("MSBB BM44").ensemblID)
    ad_msbb_mb22 = list(get_alzheimers_deg("MSBB BM22").ensemblID)
    ad_msbb_mb10 = list(get_alzheimers_deg("MSBB BM10").ensemblID)
    return {"BS_Caudate_SZ": bs_caudate_degs, "BS_DLPFC_SZ": bs_dlpfc_degs,
            "BS_Hippocampus_SZ": bs_hippo_degs, "CMC_DLPFC_SZ": cmc_dlpfc_degs,
            "BS_sACC_BD": bs_sacc_degs, "BS_Amyg_BD": bs_amyg_degs,
            "PSY_SZ": psy_sz, "PSY_ASD": psy_asd, "PSY_BD": psy_bd,
            "Parkinsons": pd_degs, "MAYO_AD": ad_mayo, "ROSMAP_AD": ad_rosmap,
            "MSBB_MD36_AD": ad_msbb_mb36, "MSBB_MD44_AD": ad_msbb_mb44,
            "MSBB_MD22_AD": ad_msbb_mb22, "MSBB_MD10_AD": ad_msbb_mb10}


def get_public_twas():
    # TWAS
    ## BrainSeq SZ
    bs_caudate_twas = list(get_bs_twas()[(get_bs_twas()["Caudate_FDR"] < 0.05)].Geneid)
    bs_dlpfc_twas = list(get_bs_twas()[(get_bs_twas()["DLPFC_FDR"] < 0.05)].Geneid)
    bs_hippo_twas = list(get_bs_twas()[(get_bs_twas()["HIPPO_FDR"] < 0.05)].Geneid)
    ## LIBD Bipolar
    bs_sacc_twas = list(get_libd_bd_twas()[(get_libd_bd_twas()["Region"] == "sACC") &
                                          (get_libd_bd_twas()["TWAS.fdrP"]<0.05)].ID)
    bs_amyg_twas = list(get_libd_bd_twas()[(get_libd_bd_twas()["Region"] == "Amygdala") &
                                          (get_libd_bd_twas()["TWAS.fdrP"]<0.05)].ID)
    ## PsychENCODE
    psy_asd_twas = list(get_gandal_twas(config["twas_asd_file"]).GeneID)
    psy_sz_twas = list(get_gandal_twas(config["twas_sz_file"]).GeneID)
    psy_bd_twas = list(get_gandal_twas(config["twas_bd_file"]).ID)
    ## Alzheimer's
    ad_twas = list(get_ad_twas().ensemblID)
    return {"BS_Caudate_SZ": bs_caudate_twas, "BS_DLPFC_SZ": bs_dlpfc_twas,
            "BS_Hippocampus_SZ": bs_hippo_twas, "PSY_SZ": psy_sz_twas,
            "PSY_ASD": psy_asd_twas, "PSY_BD": psy_bd_twas,
            "BS_sACC_BD": bs_sacc_twas, "BS_Amyg_BD": bs_amyg_twas,
            "Alzheimers": ad_twas}


def enrichment_loop(comp_list, comp_dict, label):
    dir_dict = {"all": "All", "up": "Decreased in AA",
                "down": "Increased in AA"}
    dt = pd.DataFrame()
    for direction in ["all", "up", "down"]:
        ## Load ancestry-related DEGs
        ancestry_dict = load_ancestry_de(direction)
        ## DEGs
        or_lt = []; pval_lt = []; tissue_lt = [];
        comp_lt = []; cil_lt = []; ciu_lt = []
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            for comp in comp_list:
                cil, ciu, oddratio, pvals = \
                    logreg_enrichment(ancestry_dict[tissue], comp_dict[comp])
                or_lt.append(oddratio); pval_lt.append(pvals);
                tissue_lt.append(tissue); comp_lt.append(comp);
                cil_lt.append(cil); ciu_lt.append(ciu);
        fdr = multipletests(pval_lt, method='fdr_bh')[1]
        dtx = pd.DataFrame({"Tissue": tissue_lt, "Dataset": comp_lt,
                            "OR": or_lt, "CI_Lower": cil_lt, "CI_Upper": ciu_lt,
                            "P-value": pval_lt, "FDR": fdr,
                            "Direction":dir_dict[direction], "Method":label})
        dt = pd.concat([dt, dtx], axis=0)
    return dt


def main():
    ## Load public DEGs
    degs_dict = get_public_degs()
    ## Load public TWAS
    twas_dict = get_public_twas()
    ## Labels
    degs_list = ["BS_Caudate_SZ", "BS_DLPFC_SZ", "BS_Hippocampus_SZ",
                 "BS_sACC_BD", "BS_Amyg_BD", "CMC_DLPFC_SZ", "PSY_SZ",
                 "PSY_ASD", "PSY_BD", "Parkinsons", "MAYO_AD", "ROSMAP_AD",
                 "MSBB_MD36_AD", "MSBB_MD44_AD", "MSBB_MD22_AD", "MSBB_MD10_AD"]
    twas_list = ["BS_Caudate_SZ", "BS_DLPFC_SZ", "BS_Hippocampus_SZ",
                 "BS_sACC_BD", "BS_Amyg_BD", "PSY_SZ", "PSY_ASD", "PSY_BD",
                 "Alzheimers"]
    ## Direction dictionary
    dt_degs = enrichment_loop(degs_list, degs_dict, "DEG")
    dt_twas = enrichment_loop(twas_list, twas_dict, "TWAS")
    ## Save enrichment
    dt_degs.to_csv("clincial_phenotypes_enrichment_analysis_DEGs_tx.tsv",
                   sep='\t', index=False)
    dt_twas.to_csv("clincial_phenotypes_enrichment_analysis_TWAS_tx.tsv",
                   sep='\t', index=False)
    session_info.show()


if __name__ == '__main__':
    main()
