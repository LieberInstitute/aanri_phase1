#!/bin/bash
#$ -cwd
#$ -l gpu,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N transfer_label
#$ -o ./transfer.log
#$ -e ./transfer.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

echo "**** Run tensorQTL ****"
export CUDA_VISIBLE_DEVICES=2

python3 ../_h/04.transfer_labels.py

echo "**** Job ends ****"
date
