#!/bin/bash

show_help() {
    cat <<EOF
    Usage: ${0##*/} [-h] [-t tissue user] [-f feature]
    This script extracts genotypes after filtering them.
 
    -h|--help            display this help text
    -t|--tissue  NAME    basespace user [ex: caudate]
    -f|--feature NAME    run name from basespace [ex: genes]
    -o|--outdir  NAME    output directory [default: gwas]
EOF
}

set -o errexit # stop execution on error

# Reset all variables that might be set
TISSUE=
FEATURE=
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
        -f|--feature)
            if [ -n "$2" ]; then
                FEATURE=$2
                shift
            else
                printf 'ERROR: "--feature" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --feature=?*)
            FEATURE=${1#*=} # Delete everything up to "=" and assign the remainder
            ;;
        --feature=)
            printf 'ERROR: "--feature" requires a non-empty option argument.\n' >&2
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
GENOTYPE_DIR="${BASE_LOC}/eqtl_analysis/${TISSUE}/${FEATURE}/plink_format/_m/"

mkdir -p ${TISSUE}/${FEATURE}/$OUTPUT

awk -F' ' -v OFS='\t' '{print $2,$1,$1,$2}' \
    $GENOTYPE_DIR/genotypes.fam > ${TISSUE}/${FEATURE}/$OUTPUT/id_map

for CHROM in `seq 1 22` X; do
    plink2 --bfile $GENOTYPE_DIR/genotypes \
	   --maf 0.01 --hwe 1e-5 --chr $CHROM --make-pgen \
	   --out ${TISSUE}/${FEATURE}/$OUTPUT/chr$CHROM
    
    plink2 --pfile ${TISSUE}/${FEATURE}/$OUTPUT/chr$CHROM \
	   --make-pgen --update-ids ${TISSUE}/${FEATURE}/$OUTPUT/id_map \
	   --out ${TISSUE}/${FEATURE}/$OUTPUT/chr${CHROM}_brnum
done
