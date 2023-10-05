#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N hla_enrichment
#$ -o ./summary.log
#$ -e ./summary.log

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
python3 ../_h/02_enrichment.py

echo "**** Job ends ****"
date
