#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N prepare_mashr
#$ -o ./summary.out
#$ -e ./summary.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

## Job command
FEATURE="junctions"

echo "**** Run mashr prep ****"
python ../_h/extract_top_eqtls.py --feature $FEATURE

echo "**** Job ends ****"
date
