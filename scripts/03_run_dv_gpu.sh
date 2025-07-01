#!/bin/bash

# input parameters
REFERENCE="2.no_alt_GRCh38"
CALLER="2.DV"
PREV=""
NEXT="1.filtered"
SAMPLE=$1	# ex) "raw_PGPC-02"
SOURCE=$2	# ex) "2.buccal"
GPU_NO=$3       # ex) 1
MAX_CPU=$4      # ex) 8

# base directories
DATA_DIR="/home/zunuan/anaconda3/salfilter/data"
REF_DIR="/home/zunuan/anaconda3/salfilter/reference/${REFERENCE}/WholeGenomeFasta"

BAM_DIR="${DATA_DIR}/0.raw_data/2.picard/${SAMPLE}/${SOURCE}"
INPUT_DIR="${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${SAMPLE}/${SOURCE}"
OUTPUT_DIR="${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}"

BAM=$(ls $BAM_DIR | egrep ".bam$")
ACCESSION="${BAM%.bam}"

BIN_VERSION="1.6.1"
DV_LOG_DIR="${DATA_DIR}/${REFERENCE}/${CALLER}/DV_log"


#unset DOCKER_HOST
#docker login $CI_REGISTRY -u $CI_REGISTRY_USER -p $CI_REGISTRY_PASSWORD

docker run -it --rm --device nvidia.com/gpu=${GPU_NO} --security-opt=label=disable \
  --rm --security-opt=label=disable \
  -v "${BAM_DIR}":"/input" \
  -v "${OUTPUT_DIR}":"/output" \
  -v "${REF_DIR}":"/reference" \
  -v "${DV_LOG_DIR}":"/log" \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/reference/genome.fa \
  --reads=/input/${BAM} \
  --output_vcf=/output/${ACCESSION}.vcf.gz \
  --intermediate_results_dir ${DATA_DIR}/tmp \
  --num_shards=${MAX_CPU} \
  --logging_dir=/log
