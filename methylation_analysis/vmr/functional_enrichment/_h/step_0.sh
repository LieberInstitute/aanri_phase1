#!/bin/bash

LOCAL="/dcs04/lieber/statsgen/shizhong/AANRI/VMR3/99"

cp -v $LOCAL/caud/aa2/local_ancestry/out/dmr_great_fdr0.05.csv \
   caudate_dmr_rgreat.csv

cp -v $LOCAL/dlpfc/aa2/local_ancestry/out/dmr_great_fdr0.05.csv \
   dlpfc_dmr_rgreat.csv

cp -v $LOCAL/hippo/aa2/local_ancestry/out/dmr_great_fdr0.05.csv \
   hippocampus_dmr_rgreat.csv
