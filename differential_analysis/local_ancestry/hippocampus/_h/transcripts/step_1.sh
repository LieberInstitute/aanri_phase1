#$ -cwd
#$ -R y
#$ -t 1-500:1
#$ -tc 25
#$ -l mem_free=10.0G,h_vmem=10G,h_fsize=50G
#$ -N 'local_tx_hippocampus'
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
module load pandoc

module list

echo "**** Run sensitivity analysis for DE ****"
FEATURE="transcripts"
mkdir -p $FEATURE

Rscript  ../_h/01_differential_expression.R \
	 --feature $FEATURE --sge_id $SGE_TASK_ID --chunks 500

echo "**** Job ends ****"
date -Is
