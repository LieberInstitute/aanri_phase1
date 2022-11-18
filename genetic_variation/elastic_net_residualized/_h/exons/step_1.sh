#!/bin/bash
#$ -cwd
#$ -t 1-71
#$ -tc 5
#$ -l caracol,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N predExpr_gpu_caudate_exons
#$ -o ./caudate_summary.log
#$ -e ./caudate_summary.log
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
TISSUE="caudate"; FEATURE="exons"; NSPLIT=5
START=(`seq 0 1000 70000`)
END=(`seq 1000 1000 70000` 70872)

echo "**** Run tensorQTL ****"
python ../_h/predict_expression_iter.py \
       --tissue $TISSUE --feature $FEATURE --nsplit $NSPLIT \
       --start ${START[${SGE_TASK_ID}-1]} \
       --end ${END[${SGE_TASK_ID}-1]}

echo "**** Job ends ****"
date
