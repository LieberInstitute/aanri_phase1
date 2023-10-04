#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=10G
#$ -N annotate_vmr_dlpfc
#$ -o ./dlpfc_aaOnly.log
#$ -e ./dlpfc_aaOnly.log

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
TISSUE="dlpfc"

echo "**** Annotate VMRs ****"
Rscript ../_h/vmr_annotation.R --tissue $TISSUE

echo "**** Job ends ****"
date
