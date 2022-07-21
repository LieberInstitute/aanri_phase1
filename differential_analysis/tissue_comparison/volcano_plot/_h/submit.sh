#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N volcano_plot
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
module load conda_R
module list

## Edit with your job command
echo "**** Run volcano plot ****"
Rscript ../_h/generate_volcano.R

echo "**** Job ends ****"
date
