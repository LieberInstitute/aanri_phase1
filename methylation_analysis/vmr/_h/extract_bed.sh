#!/bin/bash

## AA only
DE="../../../differential_analysis/tissue_comparison/summary_table/_m/"
mkdir -p aaOnly/
for TISSUE in Caudate DLPFC Hippocampus; do
    for FEATURE in Gene Transcript Exon Junction; do
	zcat $DE/BrainSeq_ancestry_4features_4regions.txt.gz | \
	    awk -F'\t' -v OFS='\t' -v tissue="$TISSUE" -v gene="$FEATURE" \
		'$1 == tissue && $10 == gene {print $5,$6,$7,$2}' | \
	    sort -V > aaOnly/${TISSUE,,}_de_${FEATURE,,}.bed
    done
done

## AA+EA
DE="../../../differential_analysis/internal_validation/tissue_comparison/summary_table/_m/"
mkdir -p binary/
for TISSUE in Caudate DLPFC Hippocampus; do
    for FEATURE in Gene Transcript Exon Junction; do
	zcat $DE/BrainSeq_ancestry_4features_4regions.txt.gz | \
	    awk -F'\t' -v OFS='\t' -v tissue="$TISSUE" -v gene="$FEATURE" \
		'$1 == tissue && $10 == gene {print $5,$6,$7,$2}' | \
	    sort -V > binary/${TISSUE,,}_de_${FEATURE,,}.bed
    done
done
