#!/bin/bash
#$ -cwd
#$ -t 1-86
#$ -tc 6
#$ -l caracol,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N predExpr_gpu_dentateGyrus_exons
#$ -o ./dentateGyrus_summary.log
#$ -e ./dentateGyrus_summary.log
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

export CUDA_VISIBLE_DEVICES=1

## Edit with your job command
TISSUE="dentateGyrus"; FEATURE="exons"; NSPLIT=3
NUM1=${SGE_TASK_ID}; NUM2=$((${SGE_TASK_ID} + 1))
START=${NUM1}001; END=${NUM2}000

echo "**** Run tensorQTL ****"
python ../_h/predict_expression_iter.py --tissue $TISSUE --feature $FEATURE \
       --nsplit $NSPLIT --start $START --end $END

echo "**** Job ends ****"
date
