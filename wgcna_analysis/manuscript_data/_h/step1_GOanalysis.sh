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
cp -v ../../caudate/goatools/_m/GO_analysis_module_violet.* wgcna_functional_enrichment/caudate/
cp -v ../../caudate/goatools/_m/module_annotated.csv wgcna_functional_enrichment/caudate/
# Dentate Gyrus
cp -v ../../dentateGyrus/goatools/_m/GO_analysis_module_black.* wgcna_functional_enrichment/dentateGyrus/
cp -v ../../dentateGyrus/goatools/_m/GO_analysis_module_darkgrey.* wgcna_functional_enrichment/dentateGyrus/
cp -v ../../dentateGyrus/goatools/_m/GO_analysis_module_darkolivegreen.* wgcna_functional_enrichment/dentateGyrus/
cp -v ../../dentateGyrus/goatools/_m/GO_analysis_module_lightcyan.* wgcna_functional_enrichment/dentateGyrus/
cp -v ../../dentateGyrus/goatools/_m/GO_analysis_module_magenta.* wgcna_functional_enrichment/dentateGyrus/
cp -v ../../dentateGyrus/goatools/_m/GO_analysis_module_skyblue.* wgcna_functional_enrichment/dentateGyrus/
cp -v ../../dentateGyrus/goatools/_m/module_annotated.csv wgcna_functional_enrichment/dentateGyrus/
# DLPFC
cp -v ../../dlpfc/goatools/_m/GO_analysis_module_black.* wgcna_functional_enrichment/dlpfc/
cp -v ../../dlpfc/goatools/_m/GO_analysis_module_cyan.* wgcna_functional_enrichment/dlpfc/
cp -v ../../dlpfc/goatools/_m/GO_analysis_module_darkred.* wgcna_functional_enrichment/dlpfc/
cp -v ../../dlpfc/goatools/_m/GO_analysis_module_darkturquoise.* wgcna_functional_enrichment/dlpfc/
cp -v ../../dlpfc/goatools/_m/module_annotated.csv wgcna_functional_enrichment/dlpfc/

# Hippocampus
cp -v ../../hippocampus/goatools/_m/GO_analysis_module_brown.* wgcna_functional_enrichment/hippocampus/
cp -v ../../hippocampus/goatools/_m/module_annotated.csv wgcna_functional_enrichment/hippocampus/

## Tar zip
tar -czvf WGCNA_functional_enrichment_analysis_ancestry.tar.gz wgcna_functional_enrichment
rm -rf wgcna_functional_enrichment/

echo "**** Job ends ****"
date -Is
