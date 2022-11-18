#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N extract_methylation_hippocampus
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
module load conda_R
module list

echo "**** Run methylation extraction ****"
Rscript ../_h/extract_methyl.R

echo "**** Job ends ****"
date
