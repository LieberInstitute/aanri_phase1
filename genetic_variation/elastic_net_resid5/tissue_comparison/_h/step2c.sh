#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=15.0G,h_vmem=15G,h_fsize=50G
#$ -N run_mash
#$ -o ./mash.out
#$ -e ./mash.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load R/3.6.1
module load pandoc
module list

## Edit with your job command
FEATURE="exons"

echo "**** Run mash model ****"
Rscript ../_h/mashr_script.R --feature $FEATURE

echo "**** Job ends ****"
date
