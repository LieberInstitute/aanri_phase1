#!/bin/bash
#$ -cwd
#$ -l h_fsize=10G
#$ -N vmr_corr_snps
#$ -o ./summary.log
#$ -e ./summary.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load conda_R/4.3.x
module load gcc/9.1.0
module load pandoc

module list

## Edit with your job command
echo "**** Annotate VMRs ****"
Rscript ../_h/01.correlation.R

echo "**** Job ends ****"
date
