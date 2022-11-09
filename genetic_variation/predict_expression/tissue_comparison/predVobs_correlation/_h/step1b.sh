#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -N 'effect_size_transcripts'
#$ -o ./summary_transcripts.out
#$ -e ./summary_transcripts.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load R/3.6.1
module list

echo "**** Run mash results summarization ****"
FEATURE="transcripts"
python ../_h/corr_beta.py --feature $FEATURE

echo "**** Job ends ****"
date -Is
