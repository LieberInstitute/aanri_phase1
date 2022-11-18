#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -N enet_gene_caudate
#$ -o ./logs/caudate_$TASK_ID.log
#$ -e ./logs/caudate_$TASK_ID.log
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
FEATURE="genes"; TISSUE="caudate"; KFOLD=5

echo "**** Run combine files ****"
declare -a failed=(12 25 55 91 92 99 197)

if [[ " ${failed[*]} " =~ " ${SGE_TASK_ID} " ]]; then
    Rscript ../_h/elastic_net.R \
	    --feature $FEATURE --tissue $TISSUE \
	    --sge_id $SGE_TASK_ID --k_fold $KFOLD
fi

echo "**** Job ends ****"
date
