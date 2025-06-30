#!/bin/bash

#########################
# mamba activate dotnet #
# taskset -pc 0-9 $$    #
#########################

SAMPLES=("eHOMD_yesambi_PGPC-02" "eHOMD_yesambi_PGPC-05" "eHOMD_yesambi_PGPC-06" "eHOMD_yesambi_PGPC-50")
SOURCES=("1.saliva" "2.buccal" "3.saliva_enriched" "4.buccal_enriched")

taskset -pc 0-23 $$

# Iterate over SAMPLE and SOURCE
for SAMPLE in "${SAMPLES[@]}"; do
    for SOURCE in "${SOURCES[@]}"; do
        echo "Processing SAMPLE: ${SAMPLE}, SOURCE: ${SOURCE}"

	REF_DIR="/home/zunuan/opt/canvas/Sequence/WholeGenomeFasta"
        DATA_DIR="/home/zunuan/anaconda3/salfilter/data"

        REFERENCE="2.no_alt_GRCh38"
        CALLER="CanvasAnalysis"
        PREV="0.bam"
        NEXT="1.canvas"

        # search for vcf file using *.vcf.idx in directory
        IN_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${SAMPLE}/${SOURCE}
        IN_FILE=$(ls $IN_FILE_DIR | egrep ".bam$")
        ACCESSION="${IN_FILE%.bam}"

        OUT_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}
        # make output directory if it doesnt exist
        if [ ! -d ${OUT_FILE_DIR} ]; then
          mkdir -p ${OUT_FILE_DIR}
        fi
	
	CANVAS_DIR="/home/zunuan/opt/canvas/Canvas-1.40.0.1613+master_x64"
	POP_VCF=${REF_DIR}/dbsnp.vcf
	KMER_FA=${REF_DIR}/kmer.fa
	EXCLUDE_BED=${REF_DIR}/filter13.bed # this is recommended pipeline
	
	IN_BAM=${IN_FILE_DIR}/${ACCESSION}.bam

	# make CNV.vcf.gz
	dotnet ${CANVAS_DIR}/Canvas.dll SmallPedigree-WGS \
	       --bam=${IN_BAM} \
	       --population-b-allele-vcf ${POP_VCF} \
	       -r ${KMER_FA} \
	       -g ${REF_DIR} \
	       -f ${EXCLUDE_BED} \
	       -o ${OUT_FILE_DIR}
	tabix -p vcf ${OUT_FILE_DIR}/CNV.vcf.gz

	# make CNV_noREF.vcf.gz
	zcat ${OUT_FILE_DIR}/CNV.vcf.gz | grep -v ":REF:" > ${OUT_FILE_DIR}/CNV_noREF.vcf
	bgzip ${OUT_FILE_DIR}/CNV_noREF.vcf
	tabix -p vcf ${OUT_FILE_DIR}/CNV_noREF.vcf.gz

	# make CNV_noREF_norm.vcf.gz
	bcftools norm -m -both -o ${OUT_FILE_DIR}/CNV_noREF_norm.vcf.gz -O z ${OUT_FILE_DIR}/CNV_noREF.vcf.gz
	tabix -p vcf ${OUT_FILE_DIR}/CNV_noREF_norm.vcf.gz
    done
done
