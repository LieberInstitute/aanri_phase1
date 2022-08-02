#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N extract
#$ -o ./mash.log
#$ -e ./mash.log
#$ -m e -M jade.benjamin@libd.org

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

## Job command
echo "**** Run mashr prep ****"

Rscript ../_h/extract_mash.R

echo "**** Job ends ****"
date
