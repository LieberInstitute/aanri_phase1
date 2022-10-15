#$ -cwd
#$ -R y
#$ -N 'wgcna_functional_enrichment_tar'
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

echo "**** Run script ****"

mkdir wgcna_functional_enrichment
mkdir wgcna_functional_enrichment/caudate/
mkdir wgcna_functional_enrichment/dentateGyrus/
mkdir wgcna_functional_enrichment/dlpfc/
mkdir wgcna_functional_enrichment/hippocampus/

# Caudate
cp -v ../../../caudate/goatools_analysis/_m/module* wgcna_functional_enrichment/caudate/
cp -v ../../../caudate/goatools_analysis/_m/GO* wgcna_functional_enrichment/caudate/

# Dentate Gyrus
cp -v ../../../dentateGyrus/goatools_analysis/_m/module* wgcna_functional_enrichment/dentateGyrus/
cp -v ../../../dentateGyrus/goatools_analysis/_m/GO* wgcna_functional_enrichment/dentateGyrus/

# DLPFC
cp -v ../../../dlpfc/goatools_analysis/_m/module* wgcna_functional_enrichment/dlpfc/
cp -v ../../../dlpfc/goatools_analysis/_m/GO* wgcna_functional_enrichment/dlpfc/

# Hippocampus
cp -v ../../../hippocampus/goatools_analysis/_m/module* wgcna_functional_enrichment/hippocampus/
cp -v ../../../hippocampus/goatools_analysis/_m/GO* wgcna_functional_enrichment/hippocampus/

## Tar zip
tar -czvf WGCNA_DEG_enrichment_modules_GO_analysis.tar.gz wgcna_functional_enrichment
rm -rf wgcna_functional_enrichment/

echo "**** Job ends ****"
date -Is
