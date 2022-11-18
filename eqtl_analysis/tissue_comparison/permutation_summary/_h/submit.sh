#!/bin/bash
#$ -cwd
#$ -R y
#$ -N summarize_permutation_eqtl
#$ -o ./summary.out
#$ -e ./summary.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

## Job command
echo "**** Run summarization: permutation ****"
python ../_h/plot_N_summarize.py

echo "**** Job ends ****"
date
