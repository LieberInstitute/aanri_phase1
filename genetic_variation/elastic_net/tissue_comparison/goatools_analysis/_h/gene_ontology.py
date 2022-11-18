"""
Preform gene term enrichment analysis with GOATools.
"""
import numpy as np
import session_info
import pandas as pd
from pyhere import here
from pybiomart import Dataset
from functools import lru_cache
from collections import Counter
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
def get_metrics():
    return pd.read_csv("../../metrics_summary/_m/prediction_metrics_summary.txt.gz",
                       sep='\t', compression="gzip")


@lru_cache()
def get_background(tissue):
    new_tissue = tissue.replace(" ", "_")
    return get_metrics()[(get_metrics()["region"] == new_tissue) &
                         (get_metrics()["type"] == "Gene")].copy()


@lru_cache()
def convert2entrez(tissue):
    bg = get_background(tissue)
    bg.loc[:,"ensemblID"] = bg.feature.str.replace("\\..*", "", regex=True)
    return bg.merge(get_database(), left_on='ensemblID',
                    right_on='ensembl_gene_id')


@lru_cache()
def get_de(tissue):
    fn = here("differential_analysis/tissue_comparison/summary_table/",
              "_m/BrainSeq_ancestry_4features_4regions_allFeatures.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Tissue"] == tissue) &
              (df["Type"] == "Gene")].copy()


@lru_cache()
def get_ancestry_deg(tissue, direction):
    df = get_de(tissue)
    if direction=="all":
        return df[(df["lfsr"] < 0.05)].copy()
    elif direction=="up":
        return df[(df["lfsr"] < 0.05) & (df["posterior_mean"] > 0)].copy()
    else:
        return df[(df["lfsr"] < 0.05) & (df["posterior_mean"] < 0)].copy()


@lru_cache()
def get_predicted(tissue, r2=0.01):
    df = convert2entrez(tissue)
    return df[(df["test_score_r2_median"] > r2)].copy()


@lru_cache()
def get_overlap(tissue, direction, r2):
    de_set   = set(get_ancestry_deg(tissue, direction).Effect)
    pred_set = set(get_predicted(tissue, r2).feature)
    shared_features = de_set & pred_set
    return get_predicted(tissue, r2)\
        .set_index("feature").loc[shared_features]


def obo_annotation(tissue, alpha=0.05):
    # database annotation
    bg = convert2entrez(tissue)
    fn_obo = download_go_basic_obo()
    fn_gene2go = download_ncbi_associations() # must be gunzip to work
    obodag = GODag(fn_obo) # downloads most up-to-date
    anno_hs = Gene2GoReader(fn_gene2go, taxids=[9606])
    # get associations
    ns2assoc = anno_hs.get_ns2assc()
    for nspc, id2gos in ns2assoc.items():
        print(f"{nspc} {len(id2gos):,} annotated human genes")
    goeaobj = GOEnrichmentStudyNS(
        bg['entrezgene_id'], # List of human genes with entrez IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = alpha, # default significance cut-off
        methods = ['fdr_bh'])
    return goeaobj


def run_goea(tissue, direction, r2):
    df = get_overlap(tissue, direction, r2)
    new_tissue = tissue.lower().replace(' ', '_')
    label = f"{new_tissue}_{direction}_r2_{r2}"
    geneids_study = {z[0]:z[1] for z in zip(df['entrezgene_id'],
                                            df['external_gene_name'])}
    goeaobj = obo_annotation(tissue)
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    ctr = Counter([r.NS for r in goea_results_sig])
    print(f"Significant results[{len(goea_results_sig)}] = "\
          f"{ctr['BP']} BP + {ctr['MF']} MF + {ctr['CC']} CC.")
    goeaobj.wr_xlsx("GO_analysis_mash_%s.xlsx" % label, goea_results_sig)
    goeaobj.wr_txt("GO_analysis_mash_%s.txt" % label, goea_results_sig)


def main():
    # Run GO enrichment with GOATOOLS
    for tissue in ["Caudate", "Dentate Gyrus", "DLPFC", "Hippocampus"]:
        for direction in ["all", "up", "down"]:
            for r2 in np.append(0.01, np.arange(0.05, 0.55, 0.05)):
                run_goea(tissue, direction, round(r2, 2))
    ## Session info
    session_info.show()


if __name__ == '__main__':
    main()
