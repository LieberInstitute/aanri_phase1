#$ -cwd
#$ -R y
#$ -l mem_free=40G,h_vmem=40G,h_fsize=50G
#$ -N 'go_enrichment'
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

echo "**** Run network analysis ****"
python  ../_h/gene_ontology.py

echo "**** Job ends ****"
date -Is
