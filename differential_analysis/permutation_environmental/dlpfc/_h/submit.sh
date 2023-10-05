#$ -cwd
#$ -R y
#$ -l mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -N 'permutation_dlpfc'
#$ -o ./summary.out
#$ -e ./summary.out
#$ -t 1-10:1

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load R/3.6.1
module list

echo "**** Run sensitivity analysis for DE with permutation ****"
Rscript ../_h/differential_expression.R \
	--perm_num=$SGE_TASK_ID --n_size 62

echo "**** Job ends ****"
date -Is
