# input parameters
REFERENCE="2.no_alt_GRCh38"
CALLER="2.DV"
PREV=""
NEXT="1.filtered"
SAMPLE="HROM_yesambi_PGPC-50"
SOURCE="3.saliva_enriched"
GPU_NO=3
MAX_CPU=8

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
  -v ${PWD}:"/data" \
  google/deepvariant:"${BIN_VERSION}-gpu" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/reference/genome.fa \
  --reads=/input/${BAM} \
  --output_vcf=/input/${ACCESSION}_emit-realn.vcf.gz \
  --intermediate_results_dir ${DATA_DIR}/tmp \
  --num_shards=${MAX_CPU} \
  --regions=/input/mybed.bed \
  --make_examples_extra_args="emit_realigned_reads=true,realigner_diagnostics=/input/realigned_reads"
