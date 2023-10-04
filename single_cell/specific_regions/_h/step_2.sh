#!/bin/bash
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -N convert_h5ad
#$ -o ./convert.log
#$ -e ./convert.log

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

## Edit with your job command

echo "**** Run convert r variable to H5AD ****"
Rscript ../_h/02.convert_h5ad.R

echo "**** Job ends ****"
date
