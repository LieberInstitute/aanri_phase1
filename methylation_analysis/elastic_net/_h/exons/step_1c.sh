#!/bin/bash
#$ -cwd
#$ -t 1-84
#$ -tc 45
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N predExpr_hippocampus_exons
#$ -o ./exons_summary.log
#$ -e ./exons_summary.log
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

## Edit with your job command
TISSUE="hippocampus"; FEATURE="exons"; NSPLIT=3
START=(`seq 0 1000 83000`)
END=(`seq 1000 1000 83000` 83282)

echo "**** Run predict expression ****"
python ../_h/predict_expression_iter.py \
       --tissue $TISSUE --feature $FEATURE --nsplit $NSPLIT \
       --start ${START[${SGE_TASK_ID}-1]} \
       --end ${END[${SGE_TASK_ID}-1]}

echo "**** Job ends ****"
date
