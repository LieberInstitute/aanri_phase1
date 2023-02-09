#$ -cwd
#$ -R y
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G
#$ -N 'go_analysis'
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

module list

echo "**** Run GO analysis ****"
python  ../_h/gene_ontology.py

echo "**** Job ends ****"
date -Is
