#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N extract_mhc_genes
#$ -o ./extraction.log
#$ -e ./extraction.log

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

echo "**** Run generate supplement ****"
python3 ../_h/01_extract_mhc.py

echo "**** Job ends ****"
date
