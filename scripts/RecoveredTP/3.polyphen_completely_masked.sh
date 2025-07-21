#!/bin/bash

REFERENCE="2.no_alt_GRCh38"
CALLER="5.RecoveredTP"
PREV="1-1.isec_completely_masked"
NEXT="2-1.polyphen_completely_masked"

ORAL_SOURCE_LIST=("1.saliva" "2.buccal" "3.saliva_enriched" "4.buccal_enriched")
PGPC_LIST=("PGPC-02" "PGPC-05" "PGPC-06" "PGPC-50")
for ORAL_SOURCE in "${ORAL_SOURCE_LIST[@]}"; do
  for PGPC in "${PGPC_LIST[@]}"; do
    BASE_DIR="/home/zunuan/anaconda3/salfilter/data/2.no_alt_GRCh38/5.RecoveredTP"
    
    COMPLETELY_MASKED_DIR="${BASE_DIR}/${PREV}/HROM_yesambi_${PGPC}/${ORAL_SOURCE}"
    COMPLETELY_MASKED=$(find "${COMPLETELY_MASKED_DIR}" -maxdepth 1 -type f -name "completely_masked.vcf.gz" | head -n 1)
    
    OUT_DIR="${BASE_DIR}/${NEXT}/HROM_yesambi_${PGPC}/${ORAL_SOURCE}"
    if [ ! -d ${OUT_DIR} ]; then
      mkdir -p ${OUT_DIR}
    fi
    
    vep --cache -i ${COMPLETELY_MASKED} --plugin PolyPhen_SIFT -o ${OUT_DIR}/completely_masked.txt
  done
done
