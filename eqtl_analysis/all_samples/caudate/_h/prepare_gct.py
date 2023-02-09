"""
This script is to prepare GCT files based on GTEx pipeline. This is
for interaction analysis of AA only with genetic ancestry.
"""
import argparse
import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache

def feature_map(feature):
    return {"gene": "gene", "transcript": "tx",
            "exon": "exon", "junction": "jxn"}[feature]


def to_gct(filename, df):
    description_df = pd.DataFrame({'Description': df.index.values},
                                  index=df.index)
    dfo = pd.concat([description_df, df], axis=1)
    dfo.index.name = 'Names'
    with open(filename, "wt") as out:
        print("#1.2", file=out)
        print(df.shape[0], df.shape[1], sep="\t", file=out)
        dfo.to_csv(out, sep="\t")


@lru_cache()
def get_ancestry():
    file_name = here("input/ancestry_structure/", 
                     "structure.out_ancestry_proportion_raceDemo_compare")
    return pd.read_csv(file_name, sep='\t')


@lru_cache()
def get_pheno():
    pheno_file = here("input/phenotypes/merged/_m/merged_phenotypes.csv")
    return pd.read_csv(pheno_file)\
             .merge(get_ancestry(), left_on="BrNum", right_on="id")


@lru_cache()
def get_fam():
    fam_file = here("input/genotypes/all_samples/_m/TOPMed_LIBD_AA_EA.fam")
    return pd.read_csv(fam_file, sep="\t", header=None,
                       names=["ID","BrNum","V2","V3","V4","V5"])


@lru_cache()
def load_data(feature, tissue):
    new_feature = feature_map(feature)
    pheno_df = get_pheno()[(get_pheno()["Race"].isin(["AA", "CAUC"])) &
                           (get_pheno()["Age"] > 17) &
                           (get_pheno()["Dx"].isin(["Control"]))].copy()
    pheno_df["ids"] = pheno_df.RNum
    pheno_df.set_index("ids", inplace=True)
    tpm_df = pd.read_csv(here("input/text_files_counts/tpm/_m/",
                              f"{tissue}/{feature}/tpm.csv"), index_col=0)
    counts_df = pd.read_csv(here("input/text_files_counts/_m/",
                                 f"{tissue}/{new_feature}_counts.txt"),
                            sep="\t", index_col=0)
    samples = list(set(counts_df.columns).intersection(set(pheno_df["RNum"])))
    return pheno_df.loc[samples,:], tpm_df.loc[:,samples], counts_df.loc[:,samples]


def select_idv(pheno_df, counts_df):
    samples = list(set(pheno_df.loc[counts_df.columns,:].BrNum)\
                   .intersection(set(get_fam().BrNum)))
    new_fam = get_fam()[(get_fam()["BrNum"].isin(samples))]\
        .drop_duplicates(subset="BrNum")
    new_fam.to_csv("keepFam.txt", sep='\t', index=False, header=False)
    pheno_df.loc[:, ["RNum", "BrNum", "Eur", "Afr"]]\
            .reset_index().set_index("BrNum").loc[new_fam.BrNum]\
            .to_csv("ancestry_interaction.txt", sep='\t')
    new_pheno = pheno_df.loc[:, ["RNum", "BrNum"]]\
                        .reset_index().set_index("BrNum")\
                        .loc[new_fam.BrNum].reset_index().set_index("ids")
    return new_pheno


def main():
    ## Parameters
    parser = argparse.ArgumentParser(description="Generate GCT formatted files")
    parser.add_argument("--feature", default="gene",
                        help="Feature to analyze")
    parser.add_argument("--tissue", default="caudate",
                        help="Brain region to analyze")
    args = parser.parse_args()
    feature = args.feature; tissue = args.tissue
    ## Load data
    pheno_df, tpm_df, counts_df = load_data(feature, tissue)
    ## Subset data
    new_pheno = select_idv(pheno_df, counts_df)
    genes = list(set(counts_df.index).intersection(set(tpm_df.index)))
    ## Annotate for fastQTL and/or tensorQTL
    to_gct("counts.gct", counts_df.loc[genes,new_pheno.index])
    to_gct("tpm.gct", tpm_df.loc[genes,new_pheno.index])
    ## Save files
    new_pheno.loc[:, ["RNum", "BrNum"]]\
             .to_csv("sample_id_to_brnum.tsv", sep="\t", index=False)
    pd.DataFrame({'chr':['chr'+xx for xx in [str(x) for x in range(1,23)]+['X']]})\
      .to_csv('vcf_chr_list.txt', header=False, index=None)
    ## Session information
    session_info.show()


if __name__=='__main__':
    main()
