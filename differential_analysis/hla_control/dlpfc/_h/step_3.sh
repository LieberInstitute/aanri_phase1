#$ -cwd
#$ -R y
#$ -l mem_free=15.0G,h_vmem=15G,h_fsize=50G
#$ -N 'exon_dlpfc'
#$ -o ./summary_exon.log
#$ -e ./summary_exon.log

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

echo "**** Run brain immune controlled model ****"
Rscript ../_h/differential_expression.R --feature="exons"

echo "**** Job ends ****"
date -Is
