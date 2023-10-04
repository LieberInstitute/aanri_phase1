#$ -cwd
#$ -R y
#$ -l mem_free=5.0G,h_vmem=5G,h_fsize=50G
#$ -N 'combine_tx_dlpfc'
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
FEATURE="transcripts"

Rscript  ../_h/02_combine_expr.R --feature $FEATURE

#find $FEATURE/ -type f -name '*.main' -delete
#find $FEATURE/ -type f -name '*.vmr' -delete
find $FEATURE/ -type f -name '*.r2' -delete

echo "**** Job ends ****"
date -Is
