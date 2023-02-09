#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=10G
#$ -N process_gene
#$ -o ./summary.log
#$ -e ./summary.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load gcc/9.1.0
module load R
module list

## Edit with your job command

echo "**** Run combine files ****"
Rscript ../_h/00_data_split.R --tissue "Caudate"
Rscript ../_h/00_data_split.R --tissue "Dentate Gyrus"
Rscript ../_h/00_data_split.R --tissue "DLPFC"
Rscript ../_h/00_data_split.R --tissue "Hippocampus"

mkdir caudate/logs
mkdir dentateGyrus/logs
mkdir dlpfc/logs
mkdir hippocampus/logs

echo "**** Job ends ****"
date
