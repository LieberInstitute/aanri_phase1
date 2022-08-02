#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N combine_parquet
#$ -o ./summary.log
#$ -e ./summary.log
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load samtools
module load htslib
module list

## Edit with your job command
FEATURE="genes"

echo "**** Combine parquet files ****"
python ../_h/combine_parquet.py
bgzip -f LIBD_TOPMed_AA.nominal.txt
tabix -f LIBD_TOPMed_AA.nominal.txt.gz
#echo "**** Clean directory ****"
#rm *.parquet
echo "**** Job ends ****"
date
