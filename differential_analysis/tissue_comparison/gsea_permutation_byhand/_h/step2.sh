#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-1000:1
#$ -tc 100
#$ -l mem_free=10.0G,h_vmem=10G,h_fsize=50G
#$ -N extract_mash
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
module load R
module list

## Edit with your job command
echo "**** Run extract mash ****"
Rscript ../_h/extract_effectsizes.R --perm_num=$SGE_TASK_ID 

echo "**** Job ends ****"
date
