#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=5.0G,h_vmem=5G,h_fsize=50G
#$ -N 'allele_freq'
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
module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run allele freq diff calculation ****"
Rscript ../_h/01.cal_freq_diff.R

echo "**** Job ends ****"
date
