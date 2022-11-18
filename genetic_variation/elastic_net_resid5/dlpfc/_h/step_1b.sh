#$ -cwd
#$ -R y
#$ -l mem_free=5.0G,h_vmem=5G,h_fsize=50G
#$ -N 'predicted_expr_dlpfc'
#$ -o ./summary.out
#$ -e ./summary.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load gcc/9.1.0
module load pandoc
module load R
module list

echo "**** Run sensitivity analysis for DE ****"
Rscript  ../_h/differential_expression.R --feature "transcripts"

echo "**** Job ends ****"
date -Is
