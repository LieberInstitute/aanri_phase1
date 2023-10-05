#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N immune_enrichment
#$ -o ./enrichment.log
#$ -e ./enrichment.log

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

echo "**** Run enrichment ****"
python3 ../_h/02.enrichment_fishers.py

echo "**** Job ends ****"
date
