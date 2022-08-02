# Prepare GCT files for the gtex eqtl pipeline
"""
Author: Apuã Paquola

Edited by Kynon J Benjamin

- Inputs:
    * raw counts
    * sample table
- Outputs:
    * GCT files of counts and tpm for selected samples and genes
    * A lookup table of sample_ids and brain_ids
    * A list of chromosomes to use
"""
import argparse
import functools
import pandas as pd

@functools.lru_cache()
def load_data(feature):
    fam_df = pd.read_csv("/ceph/projects/brainseq/genotype/download/topmed/"+\
                         "imputation_filter/convert2plink/filter_maf_01/_m/"+\
                         "LIBD_Brain_TopMed.fam", sep="\t", header=None,
                         names=["BrNum", "ID", "V2", "V3", "V4", "V5"])
    pheno_df = pd.read_csv("../../../../inputs/phenotypes/_m/dg_phenotypes.csv")
    pheno_df["ids"] = pheno_df.RNum
    pheno_df.set_index("ids", inplace=True)
    tpm_df = pd.read_csv("../../../../inputs/text_files_counts/tpm/_m/"+\
                         "dentateGyrus/%s/tpm.csv" % feature, index_col=0)
    counts_df = pd.read_csv("../../../../inputs/text_files_counts/_m/dentateGyrus/"+\
                            "%s_counts.txt" % feature_map(feature),
                            sep="\t", index_col=0)
    return fam_df, pheno_df, tpm_df, counts_df


def feature_map(feature):
    return {"gene": "gene", "transcript": "tx",
            "exon": "exon", "junction": "jxn"}[feature]


def to_gct(filename, df):
    description_df = pd.DataFrame({"Description": df.index.values},
                                  index=df.index)
    dfo = pd.concat([description_df, df], axis=1)
    dfo.index.name = 'Names'
    with open(filename, "wt") as out:
        print("#1.2", file=out)
        print(df.shape[0], df.shape[1], sep="\t", file=out)
        dfo.to_csv(out, sep="\t")


def select_idv(fam_df, pheno_df, counts_df):
    samples = list(set(pheno_df.loc[counts_df.columns,:].BrNum)\
                   .intersection(set(fam_df.BrNum)))
    new_fam = fam_df[(fam_df["BrNum"].isin(samples))]\
        .drop_duplicates(subset="BrNum")
    new_fam.to_csv("keepFam.txt", sep='\t', index=False, header=False)
    new_pheno = pheno_df.loc[:, ["RNum", "BrNum"]]\
                        .reset_index().set_index("BrNum")\
                        .loc[new_fam.BrNum].reset_index().set_index("ids")
    return new_pheno


def main(feature):
    fam_df, pheno_df, tpm_df, counts_df = load_data(feature)
    new_pheno = select_idv(fam_df, pheno_df, counts_df)
    genes = list(set(counts_df.index).intersection(set(tpm_df.index)))
    to_gct("counts.gct", counts_df.loc[genes,new_pheno.index])
    to_gct("tpm.gct", tpm_df.loc[genes,new_pheno.index])
    new_pheno.loc[:, ["RNum", "BrNum"]]\
             .to_csv("sample_id_to_brnum.tsv", sep="\t", index=False)
    pd.DataFrame({'chr':['chr'+xx for xx in [str(x) for x in range(1,23)]+['X']]})\
      .to_csv('vcf_chr_list.txt', header=False, index=None)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Generate GCT formatted files")
    parser.add_argument("--feature", default="gene",
                        help="Feature to analyze")
    args = parser.parse_args()
    main(args.feature)
