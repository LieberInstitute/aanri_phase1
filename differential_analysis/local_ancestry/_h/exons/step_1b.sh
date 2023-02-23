#!/bin/bash
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -N local_exon_dentateGyrus
#$ -o ./dentateGyrus/exons/logs/summary_$TASK_ID.log
#$ -e ./dentateGyrus/exons/logs/summary_$TASK_ID.log
#$ -t 1-1000
#$ -tc 50

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
FEATURE="exons"; TISSUE="dentateGyrus"

echo "**** Run combine files ****"
Rscript ../_h/01_local_ancestry.R \
	--feature $FEATURE --tissue $TISSUE \
	--sge_id $SGE_TASK_ID

echo "**** Job ends ****"
date
