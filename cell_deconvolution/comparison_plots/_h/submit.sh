#!/bin/bash
#$ -cwd
#$ -N ct_comparison
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
module load R
module list

echo "**** Run cell prop quality control ****"
Rscript ../_h/comparison_plot.R

echo "**** Job ends ****"
date
