#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N combine_hippocampus_exon
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
module load gcc/9.1.0
module load R
module list

## Edit with your job command
TISSUE="hippocampus"; FEATURE="exons"

echo "**** Run combine files ****"
Rscript ../_h/combine_files.R  --feature $FEATURE

echo "**** Job ends ****"
date
