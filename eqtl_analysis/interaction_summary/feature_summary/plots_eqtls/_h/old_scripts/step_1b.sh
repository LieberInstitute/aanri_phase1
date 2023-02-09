#!/bin/bash
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -N gwas_gene_dentateGyrus
#$ -o ./gwas.log
#$ -e ./gwas.log
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

## Edit with your job command
FEATURE="genes"

echo "**** Run combine files ****"
bash ../_h/01_gwas.sh --tissue "dentateGyrus"

echo "**** Job ends ****"
date
