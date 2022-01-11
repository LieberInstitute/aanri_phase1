#$ -cwd
#$ -R y
#$ -l mem_free=25.0G,h_vmem=30G,h_fsize=20G
#$ -N 'partial_pred_snp'
#$ -o ./log.out
#$ -e ./log.out
#$ -m e -M jade.benjamin@libd.org

## Author: Kynon Jade Benjamin

python ../_h/combine_results.py
