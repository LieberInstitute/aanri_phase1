#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-1000:1
#$ -tc 100
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -N fit_model
#$ -o ./summary.out
#$ -e ./summary.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load R
module list

## Edit with your job command
echo "**** Run model fitting ****"
TISSUE="hippocampus"; TARGET="EA"
BASELOC="/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1"
COUNTS="${BASELOC}/input/text_files_counts/_m/${TISSUE}/gene_counts.txt"
COVS="${BASELOC}/differential_analysis/${TISSUE}/extract_for_gsea/_m/genes_covariates.tsv"
ANNOT="${BASELOC}/differential_analysis/${TISSUE}/extract_for_gsea/_m/genes_annotation.txt.gz"

Rscript ../_h/fit_model.R --counts=$COUNTS --covariates=$COVS \
	--annotation=$ANNOT --target=$TARGET --perm_num=$SGE_TASK_ID \
	--output=$TISSUE

echo "**** Job ends ****"
date
