#!/bin/bash
## Author: Kynon J Benjamin
##
## This script runs the R script in parallel by tissue and feature to
## calculate partial R2. It depends on random forest results.
##
## Input: Tissue, feature, ancestry file
##
## Output: Text files for partial R2 for each individual SNP and for most
## predictive smallest subset of SNPs using random forest dRFE

TARGET="/path/to/structure.out_ancestry_proportion_raceDemo_compare"
PHENO_FILE="/path/to/phenotypes/merged_phenotypes.csv"
ML_RESULT="/path/to/dRFE/results/"
QSV_PATH="/path/to/qSV/directory/"
SNP_DIR="/path/to/onehot/directory/"

for TISSUE in 'caudate' 'dlpfc' 'hippocampus' 'dentateGyrus'; do
    mkdir -p $TISSUE
    ls -d ../../_m/$TISSUE/ENSG*/ | while read -r FILE; do
        export geneid=`basename $FILE`
        qsub -V -N partial_$TISSUE -o $TISSUE/$geneid.log -e $TISSUE/$geneid.log -cwd \
             -b y "module load R; \
                   Rscript ../_h/partial_r2.R --tissue $TISSUE --gname $geneid \
                           --ml_dir $ML_RESULT --pheno_file $PHENO_FILE \
                           --qsv_path $QSV_PATH --snp_dir $SNP_DIR --target $TARGET"
    done
done
