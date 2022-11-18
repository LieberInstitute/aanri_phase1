#!/bin/bash
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -N enet_gene_dentateGyrus
#$ -o ./logs/dentateGyrus_$TASK_ID.log
#$ -e ./logs/dentateGyrus_$TASK_ID.log
#$ -t 10-200

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
FEATURE="genes"; TISSUE="dentateGyrus"; KFOLD=3

echo "**** Run combine files ****"
declare -a failed=(25 100 196)

if [[ " ${failed[*]} " =~ " ${SGE_TASK_ID} " ]]; then
    Rscript ../_h/elastic_net.R \
	    --feature $FEATURE --tissue $TISSUE \
	    --sge_id $SGE_TASK_ID --k_fold $KFOLD
fi

echo "**** Job ends ****"
date
