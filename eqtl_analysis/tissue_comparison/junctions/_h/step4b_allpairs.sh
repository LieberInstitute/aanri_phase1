#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-1000:1
#$ -tc 50
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N mash_allpairs_array
#$ -o ./output/allpairs_$TASK_ID.log
#$ -e ./output/allpairs_$TASK_ID.log
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

## Job command
echo "**** Run mashr prep ****"

Rscript ../_h/all_association_mash.R \
	--chunk_num $SGE_TASK_ID --output output

echo "**** Job ends ****"
date
