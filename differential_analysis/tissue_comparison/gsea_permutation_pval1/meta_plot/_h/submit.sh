#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N gsea_metaplot
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load R/3.6.1
module load pandoc

module list

## Edit with your job command
echo "**** Generate immune related meta plots ****"
Rscript ../_h/metaplot_mash.R

echo "**** Job ends ****"
date
