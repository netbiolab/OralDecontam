#!/bin/bash

SAMPLES=("eHOMD_yesambi_PGPC-02" "eHOMD_yesambi_PGPC-05" "eHOMD_yesambi_PGPC-06" "eHOMD_yesambi_PGPC-50")
SOURCES=("1.saliva" "2.buccal" "3.saliva_enriched" "4.buccal_enriched")

for SAMPLE in "${SAMPLES[@]}"; do
    for SOURCE in "${SOURCES[@]}"; do
	echo "Processing SAMPLE: ${SAMPLE}, SOURCE: ${SOURCE}"
	
	REF_DIR="/home/zunuan/opt/canvas/Sequence/WholeGenomeFasta"
	DATA_DIR="/home/zunuan/anaconda3/salfilter/data"

	REFERENCE="2.no_alt_GRCh38"
	CALLER="CanvasAnalysis"
	PREV="1.canvas"
	NEXT="2.truvari"

	IN_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${SAMPLE}/${SOURCE}
        IN_FILE=$(ls $IN_FILE_DIR | egrep "CNV_noREF_norm.vcf.gz$")
        IN_ACCESSION="${IN_FILE%.vcf.gz}"


	# OUT_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}
	MOUNT_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}
	DIR_TO_MAKE=${MOUNT_DIR}/${NEXT}/${SAMPLE}
	if [ ! -d ${DIR_TO_MAKE} ]; then
	  mkdir -p ${DIR_TO_MAKE}
	fi

	PGPC=${SAMPLE: -7} # "PGPC-02"
	if [[ "$PGPC" == "PGPC-02" || "$PGPC" == "PGPC-05" || "$PGPC" == "PGPC-06" || "PGPC" == "PGPC-50" ]]; then
	  echo "PGPC match."
	fi

	GOLD_FILE_DIR=${MOUNT_DIR}/${PREV}/raw_${PGPC}/0.blood
	GOLD_FILE=$(ls $IN_FILE_DIR | egrep "CNV_noREF_norm.vcf.gz$")
	GOLD_ACCESSION="${IN_FILE%.vcf.gz}"

	MOUNTED_IN_VCF=${PREV}/${SAMPLE}/${SOURCE}/${IN_ACCESSION}.vcf.gz
	MOUNTED_GOLD_VCF=${PREV}/raw_${PGPC}/0.blood/${GOLD_ACCESSION}.vcf.gz
	MOUNTED_OUT_DIR=${NEXT}/${SAMPLE}/${SOURCE}

	docker run \
		-v ${MOUNT_DIR}:/data \
		-it truvari bench \
		-b /data/${MOUNTED_GOLD_VCF} \
		-c /data/${MOUNTED_IN_VCF} \
		-o /data/${MOUNTED_OUT_DIR}
    done
done
