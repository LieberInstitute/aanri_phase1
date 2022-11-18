#!/bin/bash
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -q shared.q@compute-0[9]*,shared.q@compute-1*&!*111
#$ -N enet_gene_dlpfc
#$ -o ./logs/gene/dlpfc_$TASK_ID.log
#$ -e ./logs/gene/dlpfc_$TASK_ID.log
#$ -t 1-200
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
module load R
module list

## Edit with your job command
FEATURE="genes"; TISSUE="dlpfc"; KFOLD=5

echo "**** Run combine files ****"
Rscript ../_h/elastic_net.R \
	--feature $FEATURE --tissue $TISSUE \
	--sge_id $SGE_TASK_ID --k_fold $KFOLD

echo "**** Job ends ****"
date
