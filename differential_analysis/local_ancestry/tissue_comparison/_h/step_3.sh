#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10.0G,h_vmem=10G,h_fsize=50G
#$ -N run_mash
#$ -o ./mash_exon.out
#$ -e ./mash_exon.out

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
module load conda_R/4.2.x

module list

## Edit with your job command
FEATURE="exons"

echo "**** Prepare and run mash modeling ****"
Rscript ../_h/extract_effectsizes.R --feature $FEATURE

echo "**** Job ends ****"
date
