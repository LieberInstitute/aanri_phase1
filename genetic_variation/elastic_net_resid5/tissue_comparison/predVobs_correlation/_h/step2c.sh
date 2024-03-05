#!/bin/bash
#SBATCH --partition=shared,bluejay
#SBATCH --job-name=corrplot_exon
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --output=plotting_exons.log

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

module load R
module list

export LD_LIBRARY_PATH="/jhpce/shared/community/core/conda_R/4.3/R/lib64/R/lib:$LD_LIBRARY_PATH"

echo "**** Run correlation plotting ****"
FEATURE="exons"
python ../_h/plotting_corr.py --feature $FEATURE

echo "**** Job ends ****"
date -Is
