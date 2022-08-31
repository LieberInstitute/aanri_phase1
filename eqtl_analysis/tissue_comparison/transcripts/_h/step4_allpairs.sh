#!/bin/bash
#$ -cwd
#$ -R y
#$ -pe local 6
#$ -l mem_free=35G,h_vmem=35G,h_fsize=100G
#$ -N mash_allpairs
#$ -o ./allpairs_mash.log
#$ -e ./allpairs_mash.log
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
	--run_chunk --chunk_size 750 --threads 6

echo "**** Job ends ****"
date
