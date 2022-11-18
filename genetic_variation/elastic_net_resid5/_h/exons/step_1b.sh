#!/bin/bash
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=10G
#$ -N gwas_exon_dentateGyrus
#$ -o ./summary.log
#$ -e ./summary.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

## Edit with your job command
FEATURE="exons"; TISSUE="dentateGyrus"

echo "**** Run combine files ****"
##bash ../_h/01_gwas.sh --feature $FEATURE --tissue $TISSUE
mkdir $TISSUE/$FEATURE/gwas
cd $TISSUE/$FEATURE/gwas
ln -s ../../genes/gwas/* .

echo "**** Job ends ****"
date
