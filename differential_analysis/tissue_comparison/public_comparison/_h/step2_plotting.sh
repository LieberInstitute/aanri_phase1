#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'public_enrichment_plotting'
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

module load R
module list

echo "**** Run network analysis ****"
Rscript  ../_h/plot_heatmap.R

echo "**** Job ends ****"
date -Is
