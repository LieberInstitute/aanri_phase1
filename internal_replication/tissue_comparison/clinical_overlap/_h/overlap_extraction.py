# Examine enrichment in psychiatric disorders TWAS and DEGs
import functools
import numpy as np
import pandas as pd
from os import environ
from pybiomart import Dataset
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '32'
config = {
    "dlpfc_file": "/ceph/users/jbenja13/phase3_paper/phase2/extract_de/_m/" +\
    "dlpfc_diffExpr_szVctl_full.txt",
    "caud8_file": "/ceph/projects/v4_phase3_paper/analysis/" +\
    "differential_expression/_m/genes/diffExpr_szVctl_full.txt",
    "hippo_file": "/ceph/users/jbenja13/phase3_paper/phase2/extract_de/_m/" +\
    "hippo_diffExpr_szVctl_full.txt",
    'cmc_file': '/ceph/projects/v3_phase3_paper/inputs/cmc/_m/' +\
    'CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500'+\
    '_gene-adjustedSVA-differentialExpression-includeAncestry-DxSCZ-DE.tsv',
    'PW_parkinson_file': '/ceph/users/jbenja13/resources/degs/' +\
    'Nido2020_PD_pwCohort.csv',
    'alzheimers_file': '/ceph/users/jbenja13/resources/degs/' +\
    'MarquesCoelho2021_AD.csv',
    "libd_bipolar_file": "/ceph/users/jbenja13/resources/degs/" +\
    "bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation_all.csv",
    'gandal_de_file': "/ceph/users/jbenja13/psychENCODE/expression_results/" +\
    "_m/gandal2018_psychENCODE_DE_results.xlsx",
    'twas_asd_file': '/ceph/users/jbenja13/psychENCODE/_m/psychENCODE_twas_asd.csv',
    'twas_sz_file': '/ceph/users/jbenja13/psychENCODE/_m/psychENCODE_twas_sz.csv',
    'twas_bd_file': '/ceph/users/jbenja13/psychENCODE/_m/psychENCODE_twas_bd.csv',
    'twas_bs_file': '/ceph/projects/v4_phase3_paper/analysis/twas/' +\
    'public_twas_comp/_m/TWAS_gene_tissue_summary.csv',
    'twas_ad_file': '/ceph/users/jbenja13/resources/twas/Gockley2021_AD.tsv',
    'twas_bd_libd_file': '/ceph/users/jbenja13/resources/twas/' +\
    'twas_2regions_bpd_output_allTested.csv',
}

@functools.lru_cache()
def get_database():
    dataset = Dataset(name="hsapiens_gene_ensembl",
                      host="http://www.ensembl.org",
                      use_cache=True)
    db = dataset.query(attributes=["ensembl_gene_id",
                                   "external_gene_name",
                                   "entrezgene_id"],
                       use_attr_names=True).dropna(subset=['entrezgene_id'])
    return db


@functools.lru_cache()
def get_deg(fn):
    return pd.read_csv(fn, sep='\t')


@functools.lru_cache()
def get_sig_deg(fn):
    return get_deg(fn)[(get_deg(fn)["adj.P.Val"] < 0.05)]


@functools.lru_cache()
def get_bd_deg():
    df = pd.read_csv(config["libd_bipolar_file"], index_col=0)
    df["ensemblID"] = df.gencodeID.str.replace("\..*", "", regex=True)
    return df[(df["Type"] == "Gene")]


@functools.lru_cache()
def get_bd_sig_deg(label):
    return get_bd_deg()[(get_bd_deg()[label] < 0.05)]


@functools.lru_cache()
def get_gandal_deg():
    return pd.read_excel(config["gandal_de_file"], sheet_name="DGE")


@functools.lru_cache()
def get_parkinson_deg():
    return pd.read_csv(config['PW_parkinson_file'])


@functools.lru_cache()
def get_alzheimers_deg(dataset):
    df = pd.read_csv(config['alzheimers_file'], skiprows=1)\
           .drop_duplicates(subset="GENEID")
    df["ensemblID"] = df.GENEID.str.replace("\..*", "", regex=True)
    return df[(df["gene.padj"] < 0.05) & (df["dataset"] == dataset)]


@functools.lru_cache()
def get_ancestry_deg(tissue, direction):
    fn = "../../../%s/_m/genes/diffExpr_EAvsAA_FDR05.txt" % tissue
    if direction=="all":
        return pd.read_csv(fn, sep='\t', index_col=0)
    elif direction=="up":
        df = pd.read_csv(fn, sep='\t', index_col=0)
        return df[(df["t"] > 0)]
    else:
        df = pd.read_csv(fn, sep='\t', index_col=0)
        return df[(df["t"] < 0)]


@functools.lru_cache()
def get_bs_twas():
    return pd.read_csv(config["twas_bs_file"])


@functools.lru_cache()
def get_libd_bd_twas():
    return pd.read_csv(config['twas_bd_libd_file'])


@functools.lru_cache()
def get_gandal_twas(fn):
    df = pd.read_csv(fn)
    return df[(df["TWAS.Bonferroni"] < 0.05)]


@functools.lru_cache()
def get_ad_twas():
    df = pd.read_csv(config['twas_ad_file'], sep='\t', header=None,
                     names=["chr", "set", "feature", "start", "end", "na1",
                            "strand", "na2", "geneid", "gene_name"])
    df["ensemblID"] = df.geneid.str.replace("\..*", "", regex=True)
    return df


def load_ancestry_degs(direction="all"):
    caudate = set(get_ancestry_deg("caudate", direction).ensemblID)
    dg = set(get_ancestry_deg("dentateGyrus", direction).ensemblID)
    dlpfc = set(get_ancestry_deg("dlpfc", direction).ensemblID)
    hippocampus = set(get_ancestry_deg("hippocampus", direction).ensemblID)
    return {"Caudate": caudate, "Dentate Gyrus": dg,
            "DLPFC": dlpfc, "Hippocampus": hippocampus}


def get_public_degs():
    # DEGs
    ## BrainSeq SZ
    bs_caudate_degs = set(get_sig_deg(config["caud8_file"]).ensemblID)
    bs_dlpfc_degs = set(get_sig_deg(config["dlpfc_file"]).ensemblID)
    bs_hippo_degs = set(get_sig_deg(config["hippo_file"]).ensemblID)
    ## CommonMind SZ
    cmc_dlpfc_degs = set(get_sig_deg(config["cmc_file"]).genes)
    ## LIBD Bipolar
    bs_sacc_degs = set(get_bd_sig_deg("adj.P.Val_sACC").ensemblID)
    bs_amyg_degs = set(get_bd_sig_deg("adj.P.Val_Amyg").ensemblID)
    ## PsychENCODE (Gandal)
    psy_sz = set(get_gandal_deg()[(get_gandal_deg()["SCZ.fdr"] < 0.05)].ensembl_gene_id)
    psy_asd = set(get_gandal_deg()[(get_gandal_deg()["ASD.fdr"] < 0.05)].ensembl_gene_id)
    psy_bd = set(get_gandal_deg()[(get_gandal_deg()["BD.fdr"] < 0.05)].ensembl_gene_id)
    ## Parkinson's
    pd_degs = set(get_parkinson_deg()[(get_parkinson_deg()["padj"] < 0.05)].EnsemblID)
    ## Alzheimer's
    ad_mayo = set(get_alzheimers_deg("MAYO").ensemblID)
    ad_rosmap = set(get_alzheimers_deg("ROSMAP").ensemblID)
    ad_msbb_mb36 = set(get_alzheimers_deg("MSBB BM36").ensemblID)
    ad_msbb_mb44 = set(get_alzheimers_deg("MSBB BM44").ensemblID)
    ad_msbb_mb22 = set(get_alzheimers_deg("MSBB BM22").ensemblID)
    ad_msbb_mb10 = set(get_alzheimers_deg("MSBB BM10").ensemblID)
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
    bs_caudate_twas = set(get_bs_twas()[(get_bs_twas()["Caudate_FDR"] < 0.05)].Geneid)
    bs_dlpfc_twas = set(get_bs_twas()[(get_bs_twas()["DLPFC_FDR"] < 0.05)].Geneid)
    bs_hippo_twas = set(get_bs_twas()[(get_bs_twas()["HIPPO_FDR"] < 0.05)].Geneid)
    ## LIBD Bipolar
    bs_sacc_twas = set(get_libd_bd_twas()[(get_libd_bd_twas()["Region"] == "sACC") &
                                          (get_libd_bd_twas()["TWAS.fdrP"]<0.05)].ID)
    bs_amyg_twas = set(get_libd_bd_twas()[(get_libd_bd_twas()["Region"] == "Amygdala") &
                                          (get_libd_bd_twas()["TWAS.fdrP"]<0.05)].ID)
    ## PsychENCODE
    psy_asd_twas = set(get_gandal_twas(config["twas_asd_file"]).GeneID)
    psy_sz_twas = set(get_gandal_twas(config["twas_sz_file"]).GeneID)
    psy_bd_twas = set(get_gandal_twas(config["twas_bd_file"]).ID)
    ## Alzheimer's
    ad_twas = set(get_ad_twas().ensemblID)
    return {"BS_Caudate_SZ": bs_caudate_twas, "BS_DLPFC_SZ": bs_dlpfc_twas,
            "BS_Hippocampus_SZ": bs_hippo_twas, "PSY_SZ": psy_sz_twas,
            "PSY_ASD": psy_asd_twas, "PSY_BD": psy_bd_twas,
            "BS_sACC_BD": bs_sacc_twas, "BS_Amyg_BD": bs_amyg_twas,
            "Alzheimers": ad_twas}


def get_overlapping_set(a, b):
    """
    Get overlap with sets a and b in universe (u). Inputs are sets.
    """
    u = set(get_database().ensembl_gene_id)
    a = set(a); b = set(b)
    yes_a = u.intersection(a)
    yes_b = u.intersection(b)
    return yes_a.intersection(yes_b)


def extract_overlap(comp_list, comp_dict, label):
    dir_dict = {"all": "All", "up": "EA Bias", "down": "AA Bias"}
    dt = pd.DataFrame()
    for direction in ["all", "up", "down"]:
        ## Load ancestry-related DEGs
        ancestry_dict = load_ancestry_degs(direction)
        ## DEGs
        tissue_lt = []; comp_lt = []
        for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
            for comp in comp_list:
                overlap = list(get_overlapping_set(ancestry_dict[tissue],
                                                   comp_dict[comp]))
                dx = pd.DataFrame({"Tissue":tissue, "Dataset": comp,
                                   "Direction": dir_dict[direction],
                                   "Method":label, "Genes": overlap})
                dt = pd.concat([dt, dx], axis=0)
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
    dt_degs = extract_overlap(degs_list, degs_dict, "DEG")
    dt_twas = extract_overlap(twas_list, twas_dict, "TWAS")
    ## Save overlap files
    df = pd.concat([dt_degs, dt_twas], axis=0)
    df.to_csv("clinical_overlap_ancestryDEGs.txt.gz", sep='\t', index=False)


if __name__ == '__main__':
    main()