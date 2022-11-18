#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=150.0G,h_vmem=150G,h_fsize=50G
#$ -N 'summarize_eqtl'
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

module list

echo "**** Run summarize results ****"
python  ../_h/feature_summarize.py

echo "**** Job ends ****"
date -Is
