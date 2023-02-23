#$ -cwd
#$ -R y
#$ -t 1-18076:1
#$ -tc 250
#$ -l mem_free=5.0G,h_vmem=5G,h_fsize=50G
#$ -N 'local_gene_dentateGyrus'
#$ -o ./summary.out
#$ -e ./summary.out

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
FEATURE="genes"
mkdir -p $FEATURE

Rscript  ../_h/01_differential_expression.R \
	 --feature $FEATURE --sge_id $SGE_TASK_ID

echo "**** Job ends ****"
date -Is
