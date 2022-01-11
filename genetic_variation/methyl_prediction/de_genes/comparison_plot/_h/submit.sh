#$ -cwd
#$ -R y
#$ -N 'comparison_plot'
#$ -o ./log.out
#$ -e ./log.out
#$ -m e -M jade.benjamin@libd.org

## Author: Kynon Jade Benjamin

module load R
Rscript ../_h/comparison_analysis.R
