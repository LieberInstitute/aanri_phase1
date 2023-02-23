#!/bin/bash
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -N chunk_gene
#$ -o ./summary.log
#$ -e ./summary.log
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load gcc/9.1.0
module load R
module list

## Edit with your job command
FEATURE="genes"; CHUNKS=200

echo "**** Run combine files ****"
Rscript ../_h/00_data_split.R --feature $FEATURE \
	--tissue "caudate" --chunks $CHUNKS
mkdir caudate/$FEATURE/logs
Rscript ../_h/00_data_split.R --feature $FEATURE \
	--tissue "dentateGyrus" --chunks $CHUNKS
mkdir dentateGyrus/$FEATURE/logs
Rscript ../_h/00_data_split.R --feature $FEATURE \
	--tissue "dlpfc" --chunks $CHUNKS
mkdir dlpfc/$FEATURE/logs
Rscript ../_h/00_data_split.R --feature $FEATURE \
	--tissue "hippocampus" --chunks $CHUNKS
mkdir hippocampus/$FEATURE/logs

echo "**** Job ends ****"
date
