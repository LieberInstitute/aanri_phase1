#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N gene_constrain
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
module load gcc/9.1.0
module load pandoc
module load R
module list

## Edit with your job command
echo "**** Run enrichment analysis ****"
python ../_h/gene_constraint.py
python ../_h/tx_constraint.py
echo "**** Plot correlation ****"
Rscript ../_h/plot_constrain.R

echo "**** Job ends ****"
date
