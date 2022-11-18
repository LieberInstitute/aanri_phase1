#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -N 'effect_size_exons'
#$ -o ./summary_exons.out
#$ -e ./summary_exons.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load R
module list

echo "**** Run mash results summarization ****"
FEATURE="exons"
python ../_h/corr_beta.py --feature $FEATURE

echo "**** Job ends ****"
date -Is
