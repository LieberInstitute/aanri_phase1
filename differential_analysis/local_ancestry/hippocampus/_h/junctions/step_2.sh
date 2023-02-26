#$ -cwd
#$ -R y
#$ -l mem_free=15.0G,h_vmem=15G,h_fsize=50G
#$ -N 'combine_jxn_hippocampus'
#$ -o ./combine.out
#$ -e ./combine.out

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
FEATURE="junctions"
Rscript  ../_h/02_combine_expr.R --feature $FEATURE

##find $FEATURE/ -type f -name '*.ancestry' -delete

echo "**** Job ends ****"
date -Is
