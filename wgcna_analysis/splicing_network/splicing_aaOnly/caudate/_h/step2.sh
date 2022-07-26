#$ -cwd
#$ -R y
#$ -pe local 15
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N 'cc_wgcna_jxn'
#$ -o ./network.out
#$ -e ./network.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load R/4.0.3
module list

echo "**** Run network analysis ****"
Rscript  ../_h/network_analysis.R 

echo "**** Job ends ****"
date -Is
