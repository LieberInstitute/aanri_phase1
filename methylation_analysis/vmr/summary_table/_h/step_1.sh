#!/bin/bash
#$ -cwd
#$ -l h_fsize=10G
#$ -N summarize_dmrs
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
module load pandoc
module list

## Edit with your job command
echo "**** Summarize DMRs (global and local) ****"

python3 ../_h/01_results_summary.py

echo "**** Job ends ****"
date
