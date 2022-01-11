#!/bin/sh
## Author: Kynon J Benjamin
##
## This script runs random forest with dynamic recusive feature elimination.
##
## Input: Onehot encoded text file of SNPs for each features and tissue.
##
## Output: dRFEtools output for random forest regression.

PHENO_FILE="/path/to/merged_phenotypes.csv"

mkdir log_files

for TISSUE in 'caudate' 'dlpfc' 'hippocampus' 'dentateGyrus'; do
    mkdir log_files/$TISSUE
    ## Replace with location for onehot text files
    ls -d ../../_m/$TISSUE/ENSG*/ | while read -r FILE; do
        export geneid=`basename $FILE`
        qsub -V -l mem_free=12.0G,h_vmem=15G -N raffe_$geneid \
             -o log_files/$TISSUE/$geneid.log -e log_files/$TISSUE/$geneid.log -cwd \
             -b y "export OPENBLAS_NUM_THREADS=1; module load gcc; module load R; \
                  python ../_h/dRFE_rf.py --tissue $TISSUE --fn $FILE \
                         --pheno_file $PHENO_FILE"
    done
done
