#!/bin/bash
#SBATCH --partition=shared,bluejay
#SBATCH --job-name=corr_genes
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=summary_genes.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility

module list

echo "**** Run mash results summarization ****"
FEATURE="genes"
python ../_h/corr_beta.py --feature $FEATURE

echo "**** Job ends ****"
date -Is
