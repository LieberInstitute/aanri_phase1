#!/bin/bash
#$ -cwd
#$ -N corr_variables
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

echo "**** Run confounder analysis ****"
Rscript ../_h/heatmap_plots.R

echo "**** Job ends ****"
date
