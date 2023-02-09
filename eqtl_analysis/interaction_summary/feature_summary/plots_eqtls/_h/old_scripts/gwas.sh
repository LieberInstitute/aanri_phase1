#!/bin/bash

show_help() {
    cat <<EOF
    Usage: ${0##*/} [-h] [-t tissue user] [-f feature]
    This script extracts genotypes after filtering them.
 
    -h|--help            display this help text
    -t|--tissue  NAME    basespace user [ex: caudate]
    -o|--outdir  NAME    output directory [default: gwas]
EOF
}

set -o errexit # stop execution on error

# Reset all variables that might be set
TISSUE=
OUTPUT="gwas"

while :; do
    case $1 in
        -h|-\?|--help)
            show_help
            exit
            ;;
        -t|--tissue)
            if [ -n "$2" ]; then
                TISSUE=$2
                shift
            else
                printf 'ERROR: "--tissue" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --tissue=?*)
            TISSUE=${1#*=} # Delete everything up to "=" and assign the remainder
            ;;
        --tissue=)
            printf 'ERROR: "--tissue" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        -o|--outdir)
            if [ -n "$2" ]; then
                OUTPUT=$2
                shift
            else
                printf 'ERROR: "--outdir" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --outdir=?*)
            OUTPUT=${1#*=}
            ;;
        --outdir=)
            printf 'ERROR: "--outdir" requires a non-empty option argument.\n' >&2
            exit 1
            ;;
        --)
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)
            break
    esac
    shift
done

BASE_LOC="/dcs04/lieber/statsgen/jbenjami/projects/aanri_phase1/"
GENOTYPE_DIR="${BASE_LOC}/eqtl_analysis/all_samples/${TISSUE}/plink_format/_m/"

mkdir -p ${TISSUE}/$OUTPUT

for CHROM in `seq 1 22` X; do
    plink2 --bfile $GENOTYPE_DIR/genotypes \
	   --maf 0.01 --hwe 1e-5 --chr $CHROM --make-pgen \
	   --out ${TISSUE}/$OUTPUT/chr${CHROM}_brnum    
done
