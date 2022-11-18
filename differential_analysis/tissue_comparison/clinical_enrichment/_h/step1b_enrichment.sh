#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'clinical_enrichment_fishers'
#$ -o ./summary.out
#$ -e ./summary.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load gcc/9.1.0
module load pandoc
module list

echo "**** Run enrichment analysis ****"
python  ../_h/clincial_enrichment_fishers_tx.py

echo "**** Job ends ****"
date -Is
