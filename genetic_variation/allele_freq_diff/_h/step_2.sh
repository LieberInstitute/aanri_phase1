#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=50.0G,h_vmem=50G,h_fsize=50G
#$ -N 'allele_freq'
#$ -o ./calc_freq.log
#$ -e ./calc_freq.log

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

echo "**** Run allele freq diff calculation ****"
python3 ../_h/02.eSNP_allele_diff.py

echo "**** Job ends ****"
date
