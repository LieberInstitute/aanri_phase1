#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'plot_go'
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

echo "**** Run GO plotting ****"
python  ../_h/plot_go.py

echo "**** Job ends ****"
date -Is
