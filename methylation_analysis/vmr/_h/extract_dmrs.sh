#!/bin/bash

VMR="/dcs04/lieber/statsgen/shizhong/AANRI/VMR3/99"

## AA only
mkdir -p aaOnly/
cp -v $VMR/caud/aa2/out/dmr.csv aaOnly/caudate_dmr.csv
cp -v $VMR/dlpfc/aa2/out/dmr.csv aaOnly/dlpfc_dmr.csv
cp -v $VMR/hippo/aa2/out/dmr.csv aaOnly/hippocampus_dmr.csv

## Local ancestry
# mkdir -p local/
cp -v $VMR/caud/aa2/local_ancestry/out/dmr.csv \
   aaOnly/caudate_dmr_local.csv
cp -v $VMR/dlpfc/aa2/local_ancestry/out/dmr.csv \
   aaOnly/dlpfc_dmr_local.csv
cp -v $VMR/hippo/aa2/local_ancestry/out/dmr.csv \
   aaOnly/hippocampus_dmr_local.csv

# ## AA + EA
# mkdir -p binary/
# cp -v $VMR/caud/aa_ea/out/dmr.csv binary/caudate_dmr.csv
# cp -v $VMR/dlpfc/aa_ea/out/dmr.csv binary/dlpfc_dmr.csv
# cp -v $VMR/hippo/aa_ea/out/dmr.csv binary/hippocampus_dmr.csv
