#!/bin/bash

VMR="/dcs04/lieber/statsgen/shizhong/AANRI/VMR3/99/"

## AA only
cp -v $VMR/caud/aa2/local_ancestry/out/vmr_residual.txt \
   aaOnly/caudate_vmr_residual.txt
cp -v $VMR/caud/aa2/local_ancestry/out/local_ancesty.csv \
   aaOnly/caudate_local_estimate.txt
cp -v $VMR/caud/aa2/out/vmr.bed aaOnly/vmr_caudate.bed

cp -v $VMR/dlpfc/aa2/local_ancestry/out/vmr_residual.txt \
   aaOnly/dlpfc_vmr_residual.txt
cp -v $VMR/dlpfc/aa2/local_ancestry/out/local_ancesty.csv \
   aaOnly/dlpfc_local_estimate.txt
cp -v $VMR/dlpfc/aa2/out/vmr.bed aaOnly/vmr_dlpfc.bed

cp -v $VMR/hippo/aa2/local_ancestry/out/vmr_residual.txt \
   aaOnly/hippocampus_vmr_residual.txt
cp -v $VMR/hippo/aa2/local_ancestry/out/local_ancesty.csv \
   aaOnly/hippocampus_local_estimate.txt
cp -v $VMR/hippo/aa2/out/vmr.bed aaOnly/vmr_hippocampus.bed
