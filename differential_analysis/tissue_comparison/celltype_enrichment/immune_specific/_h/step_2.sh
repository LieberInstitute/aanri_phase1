#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'celltype_enrichment_fishers'
#$ -o ./enrichment.log
#$ -e ./enrichment.log

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module list

echo "**** Run enrichment analysis ****"
python3  ../_h/02.celltype_enrichment.py

echo "**** Job ends ****"
date -Is
