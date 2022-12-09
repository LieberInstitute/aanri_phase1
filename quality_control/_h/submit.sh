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
module load R/3.6.1
module list

echo "**** Run general correlation ****"
Rscript ../_h/correlate_variables.R

echo "**** Job ends ****"
date
