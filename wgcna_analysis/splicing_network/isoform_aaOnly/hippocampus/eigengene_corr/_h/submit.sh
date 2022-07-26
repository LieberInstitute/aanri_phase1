#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'hh_eigen_corr'
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

module load R/4.0.3
module list

echo "**** Run network analysis ****"
Rscript  ../_h/correlate_phenotypes.R

echo "**** Job ends ****"
date -Is
