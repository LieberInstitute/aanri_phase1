#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N disease_plotting
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
module load conda_R/4.3.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run enrichment plotting ****"
Rscript ../_h/02.plot_heatmap.R

echo "**** Job ends ****"
date
