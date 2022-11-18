#!/bin/bash
#$ -cwd
#$ -l h_fsize=10G
#$ -N combine_transcript_dentateGyrus
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
module load gcc/9.1.0
module load pandoc
module load R
module list

## Edit with your job command
FEATURE="transcripts"; TISSUE="dentateGyrus"

echo "**** Run combine files ****"
Rscript ../_h/04_combine_expr.R \
	--feature $FEATURE --tissue $TISSUE

echo "**** Job ends ****"
date
