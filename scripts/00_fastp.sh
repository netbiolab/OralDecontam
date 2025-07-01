#!/bin/bash
 
# input directories
DATA_DIR=$1
CURRENT=$2
SAMPLE=$3
SOURCE=$4
NEXT=$5

# search for fastq files in directory
FQ_DIR=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}
FQ1=$(ls $FQ_DIR | egrep "*_1.fastq")  # ex) SRR123456_1.fastq
FQ2=$(ls $FQ_DIR | egrep "*_2.fastq")  # ex) SRR123456_2.fastq
ACCESSION="${FQ1%_1.fastq}" # ex) SRR123456

# make output directory if it doesnt exist
if [ ! -d ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE} ]; then
  mkdir -p ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}
fi

# parse parameters
in_fq1=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}/${ACCESSION}_1.fastq
in_fq2=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}/${ACCESSION}_2.fastq
out_fq1=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}_trimmed_1.fastq
out_fq2=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}_trimmed_2.fastq
html_log=/home/zunuan/script/log/${ACCESSION}_trimmed.html
json_log=/home/zunuan/script/log/${ACCESSION}_trimmed.json

fastp \
	-i ${in_fq1} -I ${in_fq2} \
	-o ${out_fq1} -O ${out_fq2} \
	-q 15 \
	-w 3 \
	-l 15 \
	--detect_adapter_for_pe \
	-h ${html_log} \
	-j ${json_log}
