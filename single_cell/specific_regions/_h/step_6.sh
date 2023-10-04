#!/bin/bash
#$ -cwd
#$ -l h_fsize=50G
#$ -N post_hoc
#$ -o ./tukey_hsd.log
#$ -e ./tukey_hsd.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load pandoc

module list

## Edit with your job command

echo "**** Run cell type proportions ****"
python3 ../_h/06.post_hoc.py

echo "**** Job ends ****"
date
