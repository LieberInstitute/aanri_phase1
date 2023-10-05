#$ -cwd
#$ -R y
#$ -l mem_free=15.0G,h_vmem=15G,h_fsize=50G
#$ -N 'jxn_dlpfc'
#$ -o ./summary_jxn.log
#$ -e ./summary_jxn.log

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load R/3.6.1
module list

echo "**** Run brain immune controlled model ****"
Rscript ../_h/differential_expression.R --feature="junctions"

echo "**** Job ends ****"
date -Is
