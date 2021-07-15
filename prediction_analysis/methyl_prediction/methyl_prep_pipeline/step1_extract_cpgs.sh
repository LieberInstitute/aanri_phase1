#!/bin/bash
## Author: Kynon J Benjamin
##
## This script extracts CpGs with a 20kbp window around the start position of a
## feature. Feature annotation is present within FEATURE file.
##
## Input: Annotated features, BS object with normalized WGBS data, and
## brain region.
##
## Output: CSV text file with CpG normalized methylation.

FEATURE="../../_m/degs_annotation.txt"
METHYL_DIR="/path/to/methylation/BSobject/"

for TISSUE in 'caudate' 'dlpfc'; do
    mkdir $TISSUE
    grep $TISSUE $FEATURE | while read -r PARAM; do
        GNAME=`echo $PARAM | awk '{print $1}'`
        CHR=`echo $PARAM | awk '{print $4}'`
        P0=`echo $PARAM | awk '{print $5}' | awk '{if($1<0)print 0;else print}'`
        P1=`echo $PARAM | awk '{print $6}'`
        # echo $TISSUE $GNAME $CHR $P0 $P1
        mkdir $TISSUE/$GNAME
        qsub -V -l mem_free=12.0G,h_vmem=15G -N extractDNAm_${GNAME} \
             -o $TISSUE/${GNAME}.log -e $TISSUE/${GNAME}.log -cwd \
             -b y "module load R; \
                   Rscript ../_h/extract_methyl.R --chrom $CHR --start $P0 \
                           --end $P1 --output $TISSUE/$GNAME --tissue $TISSUE \
                           --methyl_dir $METHYL_DIR"
    done
done
