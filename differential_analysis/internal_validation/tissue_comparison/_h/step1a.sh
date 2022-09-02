#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N prepare_mash
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

## Edit with your job command
FEATURE="genes"

echo "**** Prepare files ****"
mkdir $FEATURE
python ../_h/prep_files.py --feature $FEATURE

echo "**** Job ends ****"
date
