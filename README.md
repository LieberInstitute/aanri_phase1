# Genetic and environmental contributions to ancestry differences in gene expression in the human brain

This repository contains the code for the AANRI paper examining genetic
and environmental contribution to ancestry differences in postmortem brain.

**Authors:** Kynon J.M. Benjamin<sup>1-3\*</sup>, Qiang Chen<sup>1</sup>,
Nicholas J. Eagles<sup>1</sup>, Louise A. Huuki-Myers<sup>1</sup>, 
Leonardo Collado-Torres<sup>1,4</sup>, Joshua M. Stolz<sup>1</sup>, 
Geo Pertea<sup>1</sup>, Joo Heon Shin<sup>1</sup>, 
Apuã C.M. Paquola<sup>1,2</sup>, Thomas M. Hyde<sup>1-3</sup>, 
Joel E. Kleinman<sup>1,3</sup>, Andrew E. Jaffe<sup>3,5,6</sup>, 
Shizhong Han<sup>1,3,7\*</sup>, and Daniel R. Weinberger<sup>1-3,5,7\*</sup>

**Affiliations:**  
<sup>1</sup>Lieber Institute for Brain Development, Baltimore, MD, USA  
<sup>2</sup>Department of Neurology, Johns Hopkins University School of Medicine, Baltimore, MD, USA  
<sup>3</sup>Department of Psychiatry and Behavioral Sciences, Johns Hopkins University School of Medicine, Baltimore, MD, USA  
<sup>4</sup>Center for Computational Biology, Johns Hopkins University, Baltimore, MD, USA  
<sup>5</sup>Department of Neuroscience, Johns Hopkins University School of Medicine, Baltimore, MD, USA  
<sup>6</sup>Neumora Therapeutics, Watertown, MA, USA  
<sup>7</sup>Department of Genetic Medicine, Johns Hopkins University School of Medicine, Baltimore, MD, USA

<sup>\*</sup>Corresponding authors

**Abstract:**  
Ancestral differences in genomic variation are determining factors in gene regulation; however, most gene expression studies have been limited to European ancestry samples or adjusted for ancestry to identify ancestry independent associations. We instead examined the impact of genetic ancestry on gene expression and DNA methylation (DNAm) in admixed African/Black American neurotypical individuals to untangle effects of genotype and environmental factors. Ancestry-associated differentially expressed genes (DEGs), transcripts, and gene networks, while notably not implicating neurons, are enriched for genes related to immune response and vascular tissue and explain up to 26% of heritability for ischemic stroke and 27% of heritability for Parkinson’s Disease. Ancestry-associated DEGs also show general enrichment for heritability of diverse immune-related traits, but depletion for psychiatric-related traits. The cell type enrichments and direction of effects vary by brain region. These DEGs are less evolutionarily constrained and are largely explained by genetic variations; roughly 15% are predicted by DNAm variation implicating environmental exposures. We also compared Black and White Americans, confirming most of these ancestry associated DEGs. Our results highlight how environment and genetic background affect genetic ancestry differences in gene expression in the human brain and risk for brain illness.

## License

<img src="https://licensebuttons.net/l/by-nc/3.0/88x31.png" alt width="88" height="31" scale="0">
Attribution-NonCommercial: CC BY-NC

This license lets others remix, tweak, and build upon our work non-commercially as long as they acknowledge our work.

[View License Deed](https://creativecommons.org/licenses/by-nc/4.0) | [View Legal Code](https://creativecommons.org/licenses/by-nc/4.0/legalcode)

## Citation

If you use anything in this repository please cite the following pre-print: https://doi.org/10.1101/2023.03.28.534458

## Files

### Supplementary data
Supplementary data (V2): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8410191.svg)](https://doi.org/10.5281/zenodo.8410191)

### Data availability

| Description              | Location                                                                                         |
| ------------------------ | ------------------------------------------------------------------------------------------------ |
| Caudate counts           | [https://erwinpaquolalab.libd.org/caudate_eqtl/](https://erwinpaquolalab.libd.org/caudate_eqtl/) |
| DLPFC counts (total RNA) | [https://eqtl.brainseq.org/phase2/](https://eqtl.brainseq.org/phase2/)                           |
| Hippocampus counts       | [https://eqtl.brainseq.org/phase2/](https://eqtl.brainseq.org/phase2/)                           |
| Dentate gyrus counts     | [http://research.libd.org/dg_hippo_paper/data.html](http://research.libd.org/dg_hippo_paper/data.html) |
| ------------------------ | ------------------------------------------------------------------------------------------------ |
| Caudate FASTQ            | dbGaP accession |
| DLPFC FASTQ              | Globus collection `jhpce#bsp2-dlpfc`                                                             |
| Hippocampus FASTQ        | Globus collection `jhpce#bsp2-hippo`                                                             |
| Dentate gyrus FASTQ      | SRA [SRP241159](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP241159)                 |
| ------------------------ | ------------------------------------------------------------------------------------------------ |
| Genotypes                | dbGaP accession [phs000979.v3.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000979.v3.p2) |
| ------------------------ | ------------------------------------------------------------------------------------------------ |
| Caudate DNAm (CpG)       | [file 1](https://libd-wgbs-sczd.s3.amazonaws.com/CpGassays.h5)                                   |
|                          | [file 2](https://libd-wgbs-sczd.s3.amazonaws.com/CpGse.rds)                                      |
| DLPFC DNAm   (CpG)       | [batch 1](https://libd-wgbs-sczd.s3.amazonaws.com/batch1_bs_combined_DLPFC_CpG.rda)              |
|                          | [batch 2](https://libd-wgbs-sczd.s3.amazonaws.com/batch2_bs_combined_DLPFC_CpG.rda)              |
| Hippocampus DNAm (CpG)   | [batch 1](https://libd-wgbs-sczd.s3.amazonaws.com/batch1_bs_combined_Hippocampus_CpG.rda)        |
|                          | [batch 2](https://libd-wgbs-sczd.s3.amazonaws.com/batch2_bs_combined_Hippocampus_CpG.rda)        |
| ------------------------ | ------------------------------------------------------------------------------------------------ |

The caudate DNAm uses h5 file storage method. To load this data in to `R`, please use the
`loadHDF5SummarizedExperiment` function to the parent directory where you downloaded the part 1
and part 2 files. For DLPFC and hippocampus DNAm, these can be loaded into `R` with `load`.
