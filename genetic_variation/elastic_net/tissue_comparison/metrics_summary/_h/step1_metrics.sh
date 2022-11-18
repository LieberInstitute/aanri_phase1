#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N metric_summary_r2
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

module load pandoc
module load R
module list

echo "**** Run mash model ****"
Rscript ../_h/metrics_summary.R

echo "**** Job ends ****"
date
