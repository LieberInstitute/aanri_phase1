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
    'gandal_de_file': "/ceph/users/jbenja13/psychENCODE/expression_results/" +\
    "_m/gandal2018_psychENCODE_DE_results.xlsx",
    'twas_asd_file': '/ceph/users/jbenja13/psychENCODE/_m/psychENCODE_twas_asd.csv',
    'twas_sz_file': '/ceph/users/jbenja13/psychENCODE/_m/psychENCODE_twas_sz.csv',
    'twas_bd_file': '/ceph/users/jbenja13/psychENCODE/_m/psychENCODE_twas_bd.csv',
    'twas_bs_file': '/ceph/projects/v4_phase3_paper/analysis/twas/' +\
    'public_twas_comp/_m/TWAS_gene_tissue_summary.csv'
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
def get_gandal_deg():
    return pd.read_excel(config["gandal_de_file"], sheet_name="DGE")


@functools.lru_cache()
def get_ancestry_deg(tissue):
    return pd.read_csv("../../../%s/_m/genes/diffExpr_EAvsAA_FDR05.txt" % tissue,
                       sep='\t', index_col=0)


@functools.lru_cache()
def get_bs_twas():
    return pd.read_csv(config["twas_bs_file"])


@functools.lru_cache()
def get_gandal_twas(fn):
    df = pd.read_csv(fn)
    return df[(df["TWAS.Bonferroni"] < 0.05)]


def fet(a, b):
    """
    Calculates Fisher's Exact test (fet) with sets a and b in universe u.
    Inputs are sets.
    """
    u = set(get_database().ensembl_gene_id)
    a = set(a); b = set(b)
    yes_a = u.intersection(a)
    yes_b = u.intersection(b)
    no_a = u - a
    no_b = u - b
    m = [[len(yes_a.intersection(yes_b)), len(no_a.intersection(yes_b))],
         [len(yes_a.intersection(no_b)), len(no_a.intersection(no_b))]]
    return fisher_exact(m)


def load_ancestry_degs():
    caudate = set(get_ancestry_deg("caudate").ensemblID)
    dg = set(get_ancestry_deg("dentateGyrus").ensemblID)
    dlpfc = set(get_ancestry_deg("dlpfc").ensemblID)
    hippocampus = set(get_ancestry_deg("hippocampus").ensemblID)
    return {"Caudate": caudate, "Dentate Gyrus": dg,
            "DLPFC": dlpfc, "Hippocampus": hippocampus}


def get_public_data():
    # DEGs
    ## BrainSeq SZ
    bs_caudate_degs = set(get_sig_deg(config["caud8_file"]).ensemblID)
    bs_dlpfc_degs = set(get_sig_deg(config["dlpfc_file"]).ensemblID)
    bs_hippo_degs = set(get_sig_deg(config["hippo_file"]).ensemblID)
    ## CommonMind SZ
    cmc_dlpfc_degs = set(get_sig_deg(config["cmc_file"]).genes)
    ## PsychENCODE (Gandal)
    psy_sz = set(get_gandal_deg()[(get_gandal_deg()["SCZ.fdr"] < 0.05)].ensembl_gene_id)
    psy_asd = set(get_gandal_deg()[(get_gandal_deg()["ASD.fdr"] < 0.05)].ensembl_gene_id)
    psy_bd = set(get_gandal_deg()[(get_gandal_deg()["BD.fdr"] < 0.05)].ensembl_gene_id)
    # TWAS
    ## BrainSeq SZ
    bs_caudate_twas = set(get_bs_twas()[(get_bs_twas()["Caudate_FDR"] < 0.05)].Geneid)
    bs_dlpfc_twas = set(get_bs_twas()[(get_bs_twas()["DLPFC_FDR"] < 0.05)].Geneid)
    bs_hippo_twas = set(get_bs_twas()[(get_bs_twas()["HIPPO_FDR"] < 0.05)].Geneid)
    ## PsychENCODE
    psy_asd_twas = set(get_gandal_twas(config["twas_asd_file"]).GeneID)
    psy_sz_twas = set(get_gandal_twas(config["twas_sz_file"]).GeneID)
    psy_bd_twas = set(get_gandal_twas(config["twas_bd_file"]).ID)
    return {"BS_Caudate_DEG": bs_caudate_degs, "BS_DLPFC_DEG": bs_dlpfc_degs,
            "BS_Hippocampus_DEG": bs_hippo_degs, "CMC_DLPFC_DEG": cmc_dlpfc_degs,
            "PSY_SZ_DEG": psy_sz, "PSY_ASD_DEG": psy_asd, "PSY_BD_DEG": psy_bd,
            "BS_Caudate_TWAS": bs_caudate_twas, "BS_DLPFC_TWAS": bs_dlpfc_twas,
            "BS_Hippocampus_TWAS": bs_hippo_twas, "PSY_SZ_TWAS": psy_sz_twas,
            "PSY_ASD_TWAS": psy_asd_twas, "PSY_BD_TWAS": psy_bd_twas}


def main():
    ## Load ancestry-related DEGs
    degs_dict = load_ancestry_degs()
    ## Load public DEGs
    comp_dict = get_public_data()
    ## Load public TWAS
    comp_list = ["BS_Caudate_DEG", "BS_DLPFC_DEG", "BS_Hippocampus_DEG",
                 "CMC_DLPFC_DEG", "PSY_SZ_DEG", "PSY_ASD_DEG", "PSY_BD_DEG",
                 "BS_Caudate_TWAS", "BS_DLPFC_TWAS", "BS_Hippocampus_TWAS",
                 "PSY_SZ_TWAS", "PSY_ASD_TWAS", "PSY_BD_TWAS"]
    or_lt = []; pval_lt = []; tissue_lt = []; comparison_lt = [];
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        for comp in comp_list:
            oddratio, pvals = fet(degs_dict[tissue], comp_dict[comp])
            or_lt.append(oddratio); pval_lt.append(pvals);
            tissue_lt.append(tissue); comparison_lt.append(comp)
    fdr = multipletests(pval_lt, method='fdr_bh')[1]
    dt = pd.DataFrame({"Tissue": tissue_lt, "Comparison": comparison_lt,
                       "OR": or_lt, "P-value": pval_lt, "FDR": fdr})
    dt.to_csv("clincial_phenotypes_enrichment_analysis.tsv",
              sep='\t', index=False)


if __name__ == '__main__':
    main()
