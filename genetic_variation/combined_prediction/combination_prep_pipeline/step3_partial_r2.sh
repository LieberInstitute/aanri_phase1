#!/bin/bash
## Author: Kynon J Benjamin
##
## This script runs the R script in parallel by tissue and feature to
## calculate partial R2. It depends on random forest results.
##
## Input: Tissue, feature, ancestry file
##
## Output: Text files for partial R2 for each individual CpG and for most
## predictive smallest subset of CpGs and SNPs using random forest dRFE

TARGET="/path/to/structure.out_ancestry_proportion_raceDemo_compare"
PHENO_FILE="/path/to/phenotypes/merged_phenotypes.csv"
RF_RESULT="/path/to/rf/results/"
QSV_PATH="/path/to/qSV/directory/"

for TISSUE in 'caudate' 'dlpfc'; do
    mkdir $TISSUE
    qsub -V -N raffe_$TISSUE -o $TISSUE/raffe.log -e $TISSUE/raffe.log -cwd \
         -b y "module load R; Rscript ../_h/partial_r2.R --tissue $TISSUE --rf"
    ls -d ../../_m/$TISSUE/ENSG*/ | while read -r FILE; do
        export geneid=`basename $FILE`
        qsub -V -N partial2_$geneid -o $TISSUE/$geneid.log -e $TISSUE/$geneid.log -cwd \
             -b y "module load R; \
                   Rscript ../_h/partial_r2.R --tissue $TISSUE --gname $geneid \
                           --target $TARGET --qsv_path $QSV_PATH \
                           --pheno_file $PHENO_FILE --rf_dir $RF_RESULT"
    done
done
