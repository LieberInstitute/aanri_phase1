#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'summarize_partial'
#$ -o ./summary.log
#$ -e ./summary.log

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module list

echo "**** Run summarize results ****"
python3  ../_h/partial_summarize.py

echo "**** Job ends ****"
date -Is
