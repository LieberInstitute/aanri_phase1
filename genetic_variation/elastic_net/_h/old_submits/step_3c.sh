#!/bin/bash
#$ -cwd
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

export CUDA_VISIBLE_DEVICES=2

## Edit with your job command
TISSUE="dlpfc"; FEATURE="genes"; NSPLIT=5
START=1001; END=2000

echo "**** Run tensorQTL ****"
python ../_h/predict_expression_iter.py --tissue $TISSUE --feature $FEATURE \
       --nsplit $NSPLIT --start $START --end $END

echo "**** Job ends ****"
date

## This script was edited from the gpu example here:
## https://gist.github.com/lahuuki/6fa1575bb15d80aeecdcabc5e2ad5f16
