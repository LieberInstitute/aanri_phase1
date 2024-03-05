#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=5G,h_vmem=5G,h_fsize=50G
#$ -N metric_summary_r2
#$ -o ./gene_summary.out
#$ -e ./gene_summary.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load gcc/9.1.0
module load pandoc
module load R
module list

echo "**** Run mash model ****"
Rscript ../_h/metrics_summary.R --feature "genes"

echo "**** Job ends ****"
date
