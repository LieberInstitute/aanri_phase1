#!/bin/bash
## Author: Kynon J Benjamin
##
## This script converts plink SNPs into onehot encoded values. This runs in
## parallel for speed. The helper script may need to be edited based on
## file locations.
##
## Input: Features, phenotype information, residualized expression.
##
## Output: Onehot encoded text file of SNPs for each features.

PHENO_FILE="../../../input/phenotypes/merged/_m/merged_phenotypes.csv"
FEATURE="../../_m/degs_annotation.txt"
PLINK_FILES="/path/to/plink/for/each/feature" ## Location of PLINK files for each feature

for TISSUE in 'caudate' 'dlpfc' 'hippocampus' 'dentateGyrus'; do
    mkdir $TISSUE
    EXPR_FILE="../../../differential_analysis/${TISSUE}/_m/genes/residualized_expression.tsv"
    grep $TISSUE $FEATURE | while read -r PARAM; do
        GNAME=`echo $PARAM | awk '{print $1}' `
        ii=`echo $GNAME | sed 's/\./_/' `
        if [ -f "$PLINK_FILES/$TISSUE/$GNAME.fam" ];then
            mkdir $TISSUE/$ii
            qsub -V -N onehot_${GNAME} -o $TISSUE/${ii}.log -e $TISSUE/${ii}.log -cwd \
                 -b y "python ../_h/genotype_data_extractor.py \
                       --plink_file_prefix $PLINK_FILES/$TISSUE/$GNAME \
                       --dirname $TISSUE/$ii --pheno_file $PHENO_FILE \
                       --expr_file $EXPR_FILE"
        else
            echo "$GNAME has no SNPs within the region!"
        fi
    done
done
