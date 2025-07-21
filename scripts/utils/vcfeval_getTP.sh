#!/bin/bash

# input parameters
REFERENCE="2.no_alt_GRCh38"
CALLER="4.HC_default"

#########################
# Annotation Method     #
#########################
PREV="3.bcftools_norm_annot"

NEXT="RecoverdTP"

#########################
# ONLY NAME (no "raw_") #
#########################
SAMPLE=$1     # ex) "PGPC-02"

#########################
# ONLY ORAL SAMPLES!!!  #
#########################
ORAL_SOURCE=$2    # ex) "2.buccal"

#########################
# METHOD ex) HROM, HOMD #
#########################
DECONTAM_METHOD=$3     # "HROM_noambi"

MAX_CORES=$4     # "4"

# base directories
DATA_DIR="/home/zunuan/anaconda3/salfilter/data"
REF_SDF_DIR="/home/zunuan/anaconda3/salfilter/reference/${REFERENCE}/WholeGenomeFasta/genome.SDF"

# filenames
GOLD_STANDARD_DIR="${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/raw_${SAMPLE}/0.blood"
RAW_ORAL_DIR="${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/raw_${SAMPLE}/${ORAL_SOURCE}"
DECONTAM_ORAL_DIR="${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${DECONTAM_METHOD}_${SAMPLE}/${ORAL_SOURCE}"

GOLD_STANDARD=${GOLD_STANDARD_DIR}/$(ls $GOLD_STANDARD_DIR | egrep ".vcf.gz$")
RAW_ORAL=${RAW_ORAL_DIR}/$(ls $RAW_ORAL_DIR | egrep ".vcf.gz$")
DECONTAM_ORAL=${DECONTAM_ORAL_DIR}/$(ls $DECONTAM_ORAL_DIR | egrep ".vcf.gz$")

# output directories (both temporary and permanent)
TEMP_RAW_ORAL_TP="${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${DECONTAM_METHOD}_${SAMPLE}/${ORAL_SOURCE}/temp-RAW_vs_ORAL"
TEMP_RAW_DECONTAM_TP="${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${DECONTAM_METHOD}_${SAMPLE}/${ORAL_SOURCE}/temp-RAW_vs_DECONTAM"
RECOVERED_TP="${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${DECONTAM_METHOD}_${SAMPLE}/${ORAL_SOURCE}/RECOVERED_TP"
LOST_TP="${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${DECONTAM_METHOD}_${SAMPLE}/${ORAL_SOURCE}/LOST_TP"

if [ ! -d ${TEMP_RAW_ORAL_TP} ]; then
  rtg vcfeval \
    -b ${GOLD_STANDARD} \
    -c ${RAW_ORAL} \
    -t ${REF_SDF_DIR} \
    -o ${TEMP_RAW_ORAL_TP} \
    -T ${MAX_CORES}
fi

if [ ! -d ${TEMP_RAW_DECONTAM_TP} ]; then
  rtg vcfeval \
    -b ${GOLD_STANDARD} \
    -c ${DECONTAM_ORAL} \
    -t ${REF_SDF_DIR} \
    -o ${TEMP_RAW_DECONTAM_TP} \
    -T ${MAX_CORES}
fi

if [ ! -d ${RECOVERED_TP} ]; then
  rtg vcfeval \
    -b ${TEMP_RAW_ORAL_TP}/tp.vcf.gz \
    -c ${TEMP_RAW_DECONTAM_TP}/tp.vcf.gz \
    -t ${REF_SDF_DIR} \
    -o ${RECOVERED_TP} \
    -T ${MAX_CORES}
fi

if [ ! -d ${LOST_TP} ]; then
  rtg vcfeval \
    -b ${TEMP_RAW_DECONTAM_TP}/tp.vcf.gz \
    -c ${TEMP_RAW_ORAL_TP}/tp.vcf.gz \
    -t ${REF_SDF_DIR} \
    -o ${LOST_TP} \
    -T ${MAX_CORES}
fi

# cd ${RECOVERED_TP}
# mv fp.vcf.gz recovered_TP.vcf.gz
# cd ${LOST_TP}
# mv fn.vcf.gz lost_TP.vcf.gz

# rm -R ${TEMP_RAW_ORAL_TP}
# rm -R ${TEMP_RAW_DECONTAM_TP}
