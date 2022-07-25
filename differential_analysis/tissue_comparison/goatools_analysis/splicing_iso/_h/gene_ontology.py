"""
Preform gene term enrichment analysis with GOATools.
"""

import pandas as pd
import collections as cx
from pybiomart import Dataset
from functools import lru_cache
from statsmodels.stats.multitest import multipletests

# GO analysis
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

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
def get_annotation():
    base_loc = "/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/input/text_files_counts/"
    config = {
        "genes": "%s/_m/caudate/gene_annotation.tsv" % base_loc,
        "transcripts": "%s/_m/caudate/tx_annotation.tsv" % base_loc,
        "exons": "%s/_m/caudate/exon_annotation.tsv" % base_loc,
        "junctions": "%s/_m/caudate/jxn_annotation.tsv" % base_loc,
    }
    return pd.read_csv(config["transcripts"], sep='\t')\
             .loc[:, ["names", "seqnames", "gencodeID"]]


@lru_cache()
def get_mash_degs():
    return pd.read_csv("../../../_m/transcripts/lfsr_feature_4tissues.txt.gz",
                       sep='\t')\
        .rename(columns={"Dentate.Gyrus": "Dentate Gyrus"})


@lru_cache()
def subset_tissue(tissue):
    return get_mash_degs().loc[:, ["Effect", tissue]]\
                          .rename(columns={tissue: "lfsr"})


@lru_cache()
def annotate_degs(tissue):
    return subset_tissue(tissue)\
        .merge(get_annotation(), left_on="Effect", right_on="names")\
        .drop(["names", "seqnames"], axis=1)


@lru_cache()
def annot_effect_size(tissue):
    return pd.read_csv("../../../_m/transcripts/posterior_mean_feature_4tissues.txt.gz",sep='\t')\
             .rename(columns={"Dentate.Gyrus": "Dentate Gyrus"})\
             .loc[:, ["Effect", tissue]]\
             .rename(columns={tissue: "beta"})\
             .merge(annotate_degs(tissue), on="Effect")


@lru_cache()
def get_bs_specific(tissue):
    df = annot_effect_size(tissue)
    df["ensemblID"] = df.gencodeID.str.replace("\\..*", "", regex=True)
    return df


@lru_cache()
def convert2entrez(tissue):
    bg = get_bs_specific(tissue)\
        .merge(get_database(), left_on='ensemblID', right_on='ensembl_gene_id')
    deg = bg[(bg["lfsr"] < 0.05)].copy()
    #return bg.drop_duplicates("gencodeID"), deg.drop_duplicates("gencodeID")
    return bg.drop_duplicates("Effect"), deg.drop_duplicates("Effect")


def get_upregulated(df):
    return df[(df["beta"] > 0)].copy()


def get_downregulated(df):
    return df[(df["beta"] < 0)].copy()


def obo_annotation(tissue, alpha=0.05):
    # database annotation
    bg, _ = convert2entrez(tissue)
    fn_obo = download_go_basic_obo()
    fn_gene2go = download_ncbi_associations() # must be gunzip to work
    obodag = GODag(fn_obo) # downloads most up-to-date
    anno_hs = Gene2GoReader(fn_gene2go, taxids=[9606])
    # get associations
    ns2assoc = anno_hs.get_ns2assc()
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated human genes".format(NS=nspc, N=len(id2gos)))
    goeaobj = GOEnrichmentStudyNS(
        bg['entrezgene_id'], # List of human genes with entrez IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = alpha, # default significance cut-off
        methods = ['fdr_bh'])
    return goeaobj


def run_goea(tissue, direction):
    _, df0 = convert2entrez(tissue)
    label0 = tissue.lower().replace(" ", "")
    if direction == "ALL":
        df = df0
        label = label0
    elif direction == "UP":
        df = get_upregulated(df0)
        label = "%s_%s" % (label0, direction.lower())
    else:
        df = get_downregulated(df0)
        label = "%s_%s" % (label0, direction.lower())
    geneids_study = {z[0]:z[1] for z in zip(df['entrezgene_id'],
                                            df['external_gene_name'])}
    goeaobj = obo_annotation(tissue)
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    ctr = cx.Counter([r.NS for r in goea_results_sig])
    print('Significant results[{TOTAL}] = {BP} BP + {MF} MF + {CC} CC'.format(
        TOTAL=len(goea_results_sig),
        BP=ctr['BP'],  # biological_process
        MF=ctr['MF'],  # molecular_function
        CC=ctr['CC'])) # cellular_component
    goeaobj.wr_xlsx("GO_analysis_mash_%s.xlsx" % label, goea_results_sig)
    goeaobj.wr_txt("GO_analysis_mash_%s.txt" % label, goea_results_sig)


def main():
    # Run GO enrichment with GOATOOLS
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        for direction in ["ALL", "UP", "DOWN"]:
            run_goea(tissue, direction)

if __name__ == '__main__':
    main()
