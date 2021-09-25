## Combine functional enrichment analysis
import numpy as np
import pandas as pd

def map_tissue(tissue):
    return {"caudate": "Caudate", "dentateGyrus": "Dentate Gyrus",
            "dlpfc": "DLPFC", "hippocampus": "Hippocampus"}[tissue]


def annotate_GO(tissue, fn, label):
    try:
        df = pd.read_excel(fn).sort_values('p_uncorrected')
        df["Log10"] = -np.log10(df['p_fdr_bh'])
        df["Tissue"] = map_tissue(tissue)
        df["Direction_of_Effect"] = label
        return df
    except FileNotFoundError:
        pass


def extract_GO(tissue):
    config = {
        'All': '../../../%s/goatools_analysis/_m/' % tissue +\
        'GO_analysis_allDEG.xlsx',
        'AA': '../../../%s/goatools_analysis/_m/' % tissue +\
        'GO_analysis_downregulated.xlsx',
        'EA': '../../../%s/goatools_analysis/_m/' % tissue +\
        'GO_analysis_upregulated.xlsx',
    }
    go_df = []
    for direction in ["All", "AA", "EA"]:
        go_df.append(annotate_GO(tissue, config[direction], direction))
    df = pd.concat(go_df, axis=0)
    return df


def main():
    bigdf = pd.DataFrame()
    for tissue in ["caudate", "dentateGyrus", "dlpfc", "hippocampus"]:
        bigdf = pd.concat([bigdf, extract_GO(tissue)], axis=0)
    bigdf.to_csv("functional_enrichment_analysis_ancestry.txt",
                 sep='\t', index=False)


if __name__ == '__main__':
    main()
