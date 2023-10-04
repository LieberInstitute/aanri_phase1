### This script uses scvi-tools to integrate and annotate
### the COPD and IPF.

import scvi, torch
import numpy as np
import session_info
import pandas as pd
import scanpy as sc
from pyhere import here
import matplotlib.pyplot as plt
from scvi.model.utils import mde

def load_ref(celltype):
    fn_dict = {"Astro": "astrocyte",
               "Micro": "microglia",
               "Oligo": "oligo"}
    ref_data = sc.read(f"./{fn_dict[celltype]}.single_cell.h5ad")
    ## Drop cell types with less than 50 cells
    celltypes = ref_data.obs.groupby("subtype").size()[
        (ref_data.obs.groupby("subtype").size() > 50)].reset_index()
    ref_data = ref_data[(ref_data.obs.subtype.isin(celltypes.subtype))]
    sc.pp.filter_genes(ref_data, min_counts=3)
    sc.pp.filter_cells(ref_data, min_counts=3)
    ref_data.obs["tech"] = "Song"
    return ref_data


def load_tran(celltype):
    ## Load data
    adata = sc.read(f"./{celltype.lower()}.aanri_brain_regions.h5ad")
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)
    adata.obs["tech"] = "Tran"
    return adata

    
def preprocess_data(celltype):
    ref_data = load_ref(celltype)
    idata = load_tran(celltype)
    adata = idata.concatenate(ref_data)
    del idata, ref_data
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep full dimension safe
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3",
        n_top_genes=2000, layer="counts",
        batch_key="tech", subset=True,
    )
    return adata


def integrate_datasets(celltype):
    adata = preprocess_data(celltype)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="tech")
    torch.set_float32_matmul_precision('high')
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train(); vae.save("./", overwrite=True)
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])
    return adata, vae


def transfer_annotation(celltype):
    adata, vae = integrate_datasets(celltype)
    adata.obs["celltype_scanvi"] = "Unknown"
    song_mask = adata.obs["tech"] == "Song"
    adata.obs.loc[song_mask, "celltype_scanvi"] = \
        adata.obs.subtype[song_mask].values
    print(np.unique(adata.obs["celltype_scanvi"], return_counts=True))
    torch.set_float32_matmul_precision('high')
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae, adata=adata,
        unlabeled_category="Unknown",
        labels_key="celltype_scanvi",
    )
    lvae.train(max_epochs=20, n_samples_per_label=100)
    adata.obs["C_scANVI"] = lvae.predict(adata)
    adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
    adata.obsm["X_mde_scanvi"] = mde(adata.obsm["X_scANVI"])
    adata.obs.C_scANVI = pd.Categorical(
        adata.obs.C_scANVI.values,
        categories=adata.obs.subtype.cat.categories,
    )
    return adata[adata.obs.tech == "Tran"]


def main():
    ## Run analysis
    for celltype in ["Micro", "Astro", "Oligo"]:
        adata = transfer_annotation(celltype)
        adata.write(f"./{celltype.lower()}.transfer_labels.h5ad")
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
