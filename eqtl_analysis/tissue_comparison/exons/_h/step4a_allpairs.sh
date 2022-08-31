#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N mash_allpairs_chunking
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
mkdir output
Rscript ../_h/generate_chunks.R --chunk_size 1500 --output output

echo "**** Job ends ****"
date
