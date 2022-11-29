#!/bin/bash
#$ -cwd
#$ -l h_fsize=10G
#$ -N summary_exon_dlpfc
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
module load gcc/9.1.0
module load pandoc
module load R
module list

## Edit with your job command
FEATURE="exons"; TISSUE="dlpfc"

echo "**** Run combine files ****"
Rscript ../_h/03_summary.R \
	--feature $FEATURE --tissue $TISSUE

echo "**** Job ends ****"
date