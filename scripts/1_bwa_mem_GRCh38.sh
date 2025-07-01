#!/bin/bash

# input directories
DATA_DIR=$1
CURRENT=$2
SAMPLE=$3
SOURCE=$4
NEXT=$5

# constant directories
REF=/home/zunuan/reference/human/NCBI/GATK.GRCh38/BWAIndex/genome.fa

# search for fastq files in directory
FQ_DIR=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}
FQ1=$(ls $FQ_DIR | egrep "*_1.fastq$")  # ex) SRR123456_trimmed_1.fastq
FQ2=$(ls $FQ_DIR | egrep "*_2.fastq$")  # ex) SRR123456_trimmed_2.fastq
ACCESSION="${FQ1%_1.fastq}" # ex) SRR123456_trimmed

# make output directory if it doesnt exist
if [ ! -d ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE} ]; then
  mkdir -p ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}
fi

# parse parameters
in_fq1=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}/${ACCESSION}_1.fastq
in_fq2=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}/${ACCESSION}_2.fastq
out_bam=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}_GRCh38-aln.bam

# align with bwa mem
bwa mem -t 15 $REF $in_fq1 $in_fq2 | samtools sort -o $out_bam
