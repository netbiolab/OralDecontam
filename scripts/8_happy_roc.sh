#!/bin/bash
### usage : PART OF STANDARD PIPELINE ###

REFERENCE=$1   # "2.no_alt_GRCh38"
CALLER=$2      # "1.HC"
PREV=$3        # "4-1.MAF_above_0.05"
NEXT=$4        # "5-1.happy_above_0.05"
SAMPLE=$5      # "raw_PGPC-02"
SOURCE=$6      # "1.saliva"

PARAMS=$7      # --roc-filter --ci-alpha --leftshift --decompose --engine
OUT_PREFIX=$8  # roc-filter-none_ci-0.95_yesshift_nodecomp_xcmp

REF=/home/zunuan/anaconda3/salfilter/reference/${REFERENCE}/WholeGenomeFasta/genome.fa
DATA_DIR=/home/zunuan/anaconda3/salfilter/data
SCRIPT_DIR=/home/zunuan/anaconda3/salfilter/script/

# search for vcf file using *.vcf.idx in directory
IN_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${SAMPLE}/${SOURCE}
IN_FILE=$(ls $IN_FILE_DIR | egrep ".vcf.gz$")
ACCESSION="${IN_FILE%.vcf.gz}"

OUT_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}
# make output directory if it doesnt exist
if [ ! -d ${OUT_FILE_DIR} ]; then
  mkdir -p ${OUT_FILE_DIR}
fi

PGPC=${SAMPLE: -7} # "PGPC-02"
if [[ "$PGPC" == "PGPC-02" || "$PGPC" == "PGPC-05" || "$PGPC" == "PGPC-06" || "PGPC" == "PGPC-50" ]]; then
  # do nothing
else
  echo "No PGPC match."
fi

GOLD_STANDARD_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/raw_${PGPC}/0.blood
GOLD_STANDARD=$(ls $GOLD_STANDARD_DIR | egrep ".vcf.gz$")
in_vcf=${IN_FILE_DIR}/${ACCESSION}.vcf.gz


hap.py ${GOLD_STANDARD} \
       ${in_vcf} \
       -r ${REF} \
       ${PARAMS} \
       -o ${OUT_PREVIX}
       # --roc QUAL --roc-filter LowQual
