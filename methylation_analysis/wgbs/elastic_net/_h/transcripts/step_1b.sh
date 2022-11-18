#!/bin/bash
#$ -cwd
#$ -t 1-38
#$ -tc 30
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N predExpr_hippocampus_transcripts
#$ -o ./transcripts_summary.log
#$ -e ./transcripts_summary.log
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
TISSUE="hippocampus"; FEATURE="transcripts"; NSPLIT=3
START=(`seq 0 1000 37000`)
END=(`seq 1000 1000 37000` 37907)

echo "**** Run predict expression ****"
python ../_h/predict_expression_iter.py \
       --tissue $TISSUE --feature $FEATURE --nsplit $NSPLIT \
       --start ${START[${SGE_TASK_ID}-1]} \
       --end ${END[${SGE_TASK_ID}-1]}

echo "**** Job ends ****"
date
