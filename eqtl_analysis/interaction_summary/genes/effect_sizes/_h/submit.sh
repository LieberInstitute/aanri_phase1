#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=20G,h_vmem=20G,h_fsize=50G
#$ -N 'enrichment_fishers'
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
module load R
module list

echo "**** Run enrichment analysis ****"
python ../_h/corr_beta.py

echo "**** Job ends ****"
date -Is
