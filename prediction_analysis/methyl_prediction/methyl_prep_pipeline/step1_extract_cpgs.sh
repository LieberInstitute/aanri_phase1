#!/bin/bash
## Author: Kynon J Benjamin
##
## This script extracts SNPs with a 20kbp window around the start position of a gene.
## It assumes that plink2 is within PATH.
##
## Input: Annotated features, plink format of SNPs (PFILE), chromosome information,
## and sample identification to re-label for BrNum in LIBD dataset.
##
## Output: BFILE (BED/BIM/FAM) for each feature being analyzed.


FEATURE="../../_m/degs_annotation.txt" ## Replace with random genes for random analysis
GENO="../../../input/genotypes/convert2plink/_m/LIBD_Brain_TopMed"
CHR_INFO="/dcl01/lieber/apaquola/genomes/human/gencode26/GRCh38.p10.ALL/chromInfo/_m/chromInfo.txt"
SAMPLES="/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/merged_batches_topmed/usable_genotypes/Genotype_ID_2222.csv"

awk -F',' -v OFS='\t' '{print 0,$2,$1,$2}' $SAMPLES | sed '1d' | sed 's/"//g' > updated.ids

for TISSUE in 'caudate' 'dlpfc' 'hippocampus' 'dentateGyrus'; do
    mkdir $TISSUE
    grep $TISSUE $FEATURE | while read -r PARAM; do
        # Get the gene positions +/- 20 kb
        GNAME=`echo $PARAM | awk '{print $1}'`
        CHR=`echo $PARAM | awk '{print $4}' | sed 's/chr//' -`
        P0=`echo $PARAM | awk '{print $5 - 2e4}' | awk '{if($1<0)print 0;else print}'`
        P1=`echo $PARAM | awk '{print $5 + 2e4}'`
        P1=`awk -vn=$CHR -ve=$P1 '$1=="chr"n {if(e > $2)print $2;else print e}' $CHR_INFO`
        echo $TISSUE $GNAME $CHR $P0 $P1
        qsub -V -N extractSNPs_${GNAME} -o $TISSUE/${GNAME}.log -e $TISSUE/${GNAME}.log -cwd \
             -b y "plink2 --pfile $GENO --from-bp $P0 --to-bp $P1 --chr $CHR \
                          --mind 0.1 --make-bed --update-ids updated.ids --out $TISSUE/$GNAME"
    done
done

cp $SAMPLES .
