#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -N run_mash
#$ -o ./mash_jxn.out
#$ -e ./mash_jxn.out
#$ -t 1-10:1

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
FEATURE="junctions"

echo "**** Prepare and run mash modeling ****"
Rscript ../_h/extract_effectsizes.R \
	--feature $FEATURE --perm_num=$SGE_TASK_ID

echo "**** Job ends ****"
date
