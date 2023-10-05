#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -N 'effect_size_transcripts'
#$ -o ./summary_transcripts.out
#$ -e ./summary_transcripts.out

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
FEATURE="transcripts"
python ../_h/corr_beta.py --feature $FEATURE

echo "**** Job ends ****"
date -Is
