#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=50G
#$ -N chunk_jxn
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

## Edit with your job command
FEATURE="junctions"; CHUNKS=1250

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
