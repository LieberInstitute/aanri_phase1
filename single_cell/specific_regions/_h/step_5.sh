#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=50G
#$ -N prop_transfer
#$ -o ./proportion.log
#$ -e ./proportion.log

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

echo "**** Run cell type proportions ****"
Rscript ../_h/05.cell_proportion.R

echo "**** Job ends ****"
date
