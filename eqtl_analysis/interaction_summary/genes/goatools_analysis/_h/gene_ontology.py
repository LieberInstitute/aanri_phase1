"""
Preform gene term enrichment analysis with GOATools.
"""
import session_info
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
def load_eqtl():
    fn = "../../_m/lfsr_allpairs_ancestry.txt.gz"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def subset_tissue(tissue):
    return load_eqtl().loc[:, ["gene_id", "variant_id", tissue]]\
                      .rename(columns={tissue: "lfsr"})


@lru_cache()
def get_bs_specific(tissue):
    df = subset_tissue(tissue)
    df["ensemblID"] = df.gene_id.str.replace("\\..*", "", regex=True)
    return df


@lru_cache()
def convert2entrez(tissue):
    return get_bs_specific(tissue)\
        .merge(get_database(), left_on='ensemblID', right_on='ensembl_gene_id')


@lru_cache()
def get_mash_eqtl(tissue):
    return convert2entrez(tissue)[(convert2entrez(tissue)["lfsr"] < 0.05)].copy()


def obo_annotation(tissue, alpha=0.05):
    # database annotation
    bg = convert2entrez(tissue).drop_duplicates(subset="gene_id")
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


@lru_cache()
def overlapping_features(tissue):
    deg = set(get_ancestry_deg(tissue).Effect)
    df = get_mash_eqtl(tissue); egene = set(df.gene_id)
    shared_features = deg & egene
    #print("There are %d overlapping features!" % len(shared_features))
    return df[(df["gene_id"].isin(shared_features))]\
        .drop_duplicates(subset="gene_id")


def run_goea(tissue):
    label = tissue.lower()
    df = get_mash_eqtl(tissue).drop_duplicates(subset="gene_id")
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
    # goeaobj.wr_xlsx("GO_analysis_mash_%s.xlsx" % label, goea_results_sig)
    # goeaobj.wr_txt("GO_analysis_mash_%s.txt" % label, goea_results_sig)
    return goea_results_sig


def main():
    # Run GO enrichment with GOATOOLS
    for tissue in ["Caudate", "Dentate_Gyrus", "DLPFC", "Hippocampus"]:
        run_goea(tissue)
    # Reproducibility information
    session_info.show()


if __name__ == '__main__':
    main()
