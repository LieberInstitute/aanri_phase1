#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=100G
#$ -N plotting
#$ -o ./plotting.log
#$ -e ./plotting.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run generate plot ****"
Rscript ../_h/03_plotting.R

echo "**** Job ends ****"
date
