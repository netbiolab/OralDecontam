#!/bin/bash

REFERENCE="2.no_alt_GRCh38"
CALLER="5.RecoveredTP"
PREV="0.bcftools_norm_annot"
NEXT="1-1.isec_completely_masked"

ORAL_SOURCE_LIST=("1.saliva" "2.buccal" "3.saliva_enriched" "4.buccal_enriched")
PGPC_LIST=("PGPC-02" "PGPC-05" "PGPC-06" "PGPC-50")
for ORAL_SOURCE in "${ORAL_SOURCE_LIST[@]}"; do
  for PGPC in "${PGPC_LIST[@]}"; do
    BASE_DIR="/home/zunuan/anaconda3/salfilter/data/2.no_alt_GRCh38/5.RecoveredTP"
    
    RAW_BLOOD_DIR="${BASE_DIR}/${PREV}/raw_${PGPC}/0.blood"
    RAW_ORAL_DIR="${BASE_DIR}/${PREV}/raw_${PGPC}/${ORAL_SOURCE}"
    DECONTAM_ORAL_DIR="${BASE_DIR}/${PREV}/HROM_yesambi_${PGPC}/${ORAL_SOURCE}"
    
    RAW_BLOOD=$(find "${RAW_BLOOD_DIR}" -maxdepth 1 -type f -name "SRR*.vcf.gz" | head -n 1)
    RAW_ORAL=$(find "${RAW_ORAL_DIR}" -maxdepth 1 -type f -name "SRR*.vcf.gz" | head -n 1)
    DECONTAM_ORAL=$(find "${DECONTAM_ORAL_DIR}" -maxdepth 1 -type f -name "SRR*.vcf.gz" | head -n 1)
    
    OUT_DIR="${BASE_DIR}/${NEXT}/HROM_yesambi_${PGPC}/${ORAL_SOURCE}"
    if [ ! -d ${OUT_DIR} ]; then
      mkdir -p ${OUT_DIR}
    fi
    
    # -e - : variants are not even present in the raw oral VCF
    bcftools isec \
	    -e 'FILTER="RefCall"' \
	    -e 'FILTER="RefCall"' \
	    -e - \
	    -n~110 \
	    -c all \
	    ${RAW_BLOOD} \
	    ${DECONTAM_ORAL} \
	    ${RAW_ORAL} \
	    -p ${OUT_DIR}

    NEEDED_ISEC_OUTPUT=${OUT_DIR}/0001.vcf
    COMPLETELY_MASKED_FROM_DECONTAM_ORAL=${OUT_DIR}/completely_masked.vcf.gz
    bgzip -c ${NEEDED_ISEC_OUTPUT} > ${COMPLETELY_MASKED_FROM_DECONTAM_ORAL}
    tabix -p vcf ${COMPLETELY_MASKED_FROM_DECONTAM_ORAL}
  done
done
