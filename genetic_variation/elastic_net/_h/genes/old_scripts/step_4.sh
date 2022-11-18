#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N predExpr_gpu_hippocampus_genes
#$ -o ./hippocampus_summary.log
#$ -e ./hippocampus_summary.log
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

export CUDA_VISIBLE_DEVICES=3

## Edit with your job command
TISSUE="hippocampus"; FEATURE="genes";
NSPLIT=5; END=16653

echo "**** Run tensorQTL ****"
python ../_h/predict_expression_loop.py \
       --tissue $TISSUE --feature $FEATURE \
       --nsplit $NSPLIT --end $END

echo "**** Job ends ****"
date
