#!/bin/bash
#$ -cwd
#$ -t 1-17
#$ -tc 5
#$ -l caracol,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N predExpr_gpu_dlpfc_genes
#$ -o ./dlpfc_summary.log
#$ -e ./dlpfc_summary.log
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load R/3.6.1
module list

export CUDA_VISIBLE_DEVICES=0

## Edit with your job command
TISSUE="dlpfc"; FEATURE="genes"; NSPLIT=5
START=(`seq 0 1000 16000`)
END=(`seq 1000 1000 16000` 16610)

echo "**** Run tensorQTL ****"
python ../_h/predict_expression_iter.py \
       --tissue $TISSUE --feature $FEATURE --nsplit $NSPLIT \
       --start ${START[${SGE_TASK_ID}-1]} \
       --end ${END[${SGE_TASK_ID}-1]}

echo "**** Job ends ****"
date
