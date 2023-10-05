#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -N 'effect_size_exons'
#$ -o ./summary_exons.out
#$ -e ./summary_exons.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.3.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run mash results summarization ****"
FEATURE="exons"
python ../_h/corr_beta.py --feature $FEATURE

echo "**** Job ends ****"
date -Is
