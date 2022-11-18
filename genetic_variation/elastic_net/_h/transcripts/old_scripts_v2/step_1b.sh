#!/bin/bash
#$ -cwd
#$ -t 1-14
#$ -tc 4
#$ -l caracol,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N predExpr_gpu_caudate_transcripts
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

export CUDA_VISIBLE_DEVICES=0

## Edit with your job command
TISSUE="caudate"; FEATURE="transcripts"; NSPLIT=5
START=(`seq 5000 2500 37500`); ##echo ${START[13-1]}
END=(`seq 7500 2500 37500` 37869)

echo "**** Run tensorQTL ****"
python ../_h/predict_expression_iter.py \
       --tissue $TISSUE --feature $FEATURE --nsplit $NSPLIT \
       --start ${START[${SGE_TASK_ID}-1]} \
       --end ${END[${SGE_TASK_ID}-1]}

echo "**** Job ends ****"
date
