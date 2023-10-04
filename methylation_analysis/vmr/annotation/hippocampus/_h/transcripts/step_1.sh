#$ -cwd
#$ -R y
#$ -t 1-50:1
#$ -tc 25
#$ -l mem_free=5G,h_vmem=5G,h_fsize=50G
#$ -N 'vmr_tx_hippocampus'
#$ -o ./summary_tx.out
#$ -e ./summary_tx.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
module load gcc/9.1.0
module load python
module load pandoc

module list

echo "**** Run sensitivity analysis for DE ****"
FEATURE="transcripts"

Rscript  ../_h/01_differential_expression.R \
	 --feature $FEATURE --sge_id $SGE_TASK_ID --chunks 50

echo "**** Job ends ****"
date -Is
