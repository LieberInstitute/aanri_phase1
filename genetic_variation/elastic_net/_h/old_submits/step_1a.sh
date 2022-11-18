#!/bin/bash
#$ -cwd
#$ -t 1-16590:1
#$ -tc 4
#$ -l caracol,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N predExpr_gpu_caudate_gene
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

USAGE_CUTOFF=10
NUM_GPUS=1

avail_gpus=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader | \
                 cut -d " " -f 1 | \
                 awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}')

#  Simply exit with an error if there are no GPUs left
if [[ -z $avail_gpus ]]; then
    echo "No GPUs are available."
    exit 1
fi

export CUDA_VISIBLE_DEVICES=$(echo "$avail_gpus" | head -n $NUM_GPUS | \
                                  paste -sd ",")

## Edit with your job command
TISSUE="caudate"
FEATURE="genes"
NSPLIT=5

echo "**** Run tensorQTL ****"
python ../_h/predict_expression.py --tissue $TISSUE --feature $FEATURE \
       --nsplit $NSPLIT --sge_id $SGE_TASK_ID

echo "**** Job ends ****"
date

## This script was edited from the gpu example here:
## https://gist.github.com/lahuuki/6fa1575bb15d80aeecdcabc5e2ad5f16
