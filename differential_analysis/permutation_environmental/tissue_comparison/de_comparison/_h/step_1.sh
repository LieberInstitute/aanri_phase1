#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -N 'comparison_de'
#$ -o ./comparison.log
#$ -e ./comparison.log

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module list

echo "**** Comparison of Black and Black v White results ****"
python3 ../_h/01_comparison.py

echo "**** Job ends ****"
date -Is
