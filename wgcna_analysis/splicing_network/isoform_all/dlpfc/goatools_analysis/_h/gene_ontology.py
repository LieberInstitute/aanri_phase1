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
def get_significant_modules():
    return pd.read_csv("../../eigengene_corr/_m/eigen_correlation_ancestry.tsv",
                       sep="\t")


@lru_cache()
def get_wgcna_modules():
    df = pd.read_csv("../../_m/modules.csv", index_col=0)
    modules = get_significant_modules()\
        .loc[get_significant_modules().Pvalue < 0.05, "Modules"]\
        .str.replace("ME", "")
    return df[(df["module"].isin(modules))].copy()


@lru_cache()
def get_annotation():
    base_loc = "/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/input/text_files_counts/"
    config = {
        "genes": "%s/_m/dlpfc/gene_annotation.tsv" % base_loc,
        "transcripts": "%s/_m/dlpfc/tx_annotation.tsv" % base_loc,
        "exons": "%s/_m/dlpfc/exon_annotation.tsv" % base_loc,
        "junctions": "%s/_m/dlpfc/jxn_annotation.tsv" % base_loc,
    }
    return pd.read_csv(config["transcripts"], sep='\t')\
             .loc[:, ["names", "seqnames", "gencodeID"]]


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
def get_background():
    df = pd.read_csv("../../_m/modules.csv", index_col=0)\
             .merge(get_annotation(), left_index=True, right_on="names")
    df["ensembl_gene_id"] = df.gencodeID.str.replace("\\..*", "", regex=True)
    return df.merge(get_database(), on="ensembl_gene_id")\
             .drop_duplicates(subset="names")\
             .drop("module", axis=1)


def convert2entrez(mod):
    return get_wgcna_modules()[(get_wgcna_modules().module) == mod]\
        .merge(get_background(), left_index=True, right_on="names")


def obo_annotation(alpha=0.05):
    # database annotation
    fn_obo = download_go_basic_obo()
    fn_gene2go = download_ncbi_associations() # must be gunzip to work
    obodag = GODag(fn_obo) # downloads most up-to-date
    anno_hs = Gene2GoReader(fn_gene2go, taxids=[9606])
    # get associations
    ns2assoc = anno_hs.get_ns2assc()
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated human genes".format(NS=nspc, N=len(id2gos)))
    goeaobj = GOEnrichmentStudyNS(
        get_background()['entrezgene_id'], # List of human genes with entrez IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = alpha, # default significance cut-off
        methods = ['fdr_bh'])
    return goeaobj


def run_goea(mod):
    df = convert2entrez(mod)
    geneids_study = {z[0]:z[1] for z in zip(df['entrezgene_id'],
                                            df['external_gene_name'])}
    goeaobj = obo_annotation()
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    ctr = cx.Counter([r.NS for r in goea_results_sig])
    print('Significant results[{TOTAL}] = {BP} BP + {MF} MF + {CC} CC'.format(
        TOTAL=len(goea_results_sig),
        BP=ctr['BP'],  # biological_process
        MF=ctr['MF'],  # molecular_function
        CC=ctr['CC'])) # cellular_component
    goeaobj.wr_xlsx("GO_analysis_module_%s.xlsx" % mod, goea_results_sig)
    goeaobj.wr_txt("GO_analysis_module_%s.txt" % mod, goea_results_sig)


def main():
    # Annotate WGCNA significant modules
    get_wgcna_modules().merge(get_background(), left_index=True, right_on="names")\
                       .to_csv("modules_annotated.csv", index=False)
    # Run GO enrichment with GOATOOLS for each module
    for mod in get_wgcna_modules().module.unique():
        run_goea(mod)


if __name__ == '__main__':
    main()
