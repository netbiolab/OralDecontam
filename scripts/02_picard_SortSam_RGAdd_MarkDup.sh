#!/bin/bash

# input directories
DATA_DIR=$1
CURRENT=$2
SAMPLE=$3
SOURCE=$4
NEXT=$5

# constant directories
PICARD_DIR=/home/zunuan/opt/picard
TEMP_DIR=/home/zunuan/tmp

# search for sam file in directory
BAM_DIR=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}
BAM=$(ls $BAM_DIR | egrep ".bam$")
ACCESSION="${BAM%.bam}"

# make output directory if it doesnt exist
if [ ! -d ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE} ]; then
  mkdir -p ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}
fi

# parse parameters
in_bam=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}/${ACCESSION}.bam
inter_bam_RGadded=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}_RGadded.tmp.bam
out_bam_RGadded_dedup=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}_dedup.bam
metrics=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}_dedup.metrics

/usr/bin/java -XX:ParallelGCThreads=2 -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups \
        INPUT=$in_bam \
        OUTPUT=$inter_bam_RGadded \
        SORT_ORDER=coordinate \
        RGLB='LIB1' \
        RGPL='illumina' \
        RGPU='WGS' \
        RGSM=$SAMPLE \
        CREATE_INDEX=false \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$TEMP_DIR

# @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1
# ID : Lane identifier or instrument if not applicable.
# SM : Sample identifier (ex. NA12828)
# PL : Platform used (ex. illumina, pacbio, iontorrent)
# LB : Library from which the DNA was sequenced
# PU : Platform unit identifier.

# no multithread options
/usr/bin/java -XX:ParallelGCThreads=2 -jar $PICARD_DIR/picard.jar MarkDuplicates \
        INPUT=$inter_bam_RGadded \
        OUTPUT=$out_bam_RGadded_dedup \
        METRICS_FILE=$metrics \
        CREATE_INDEX=true \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$TEMP_DIR

# remove intermediates
rm $inter_bam_RGadded
