#!/bin/bash

REFERENCE="2.no_alt_GRCh38"
CALLER="2.DV"
PREV="RecoverdTP_parsed"
NEXT="RecoverdTP_parsed_clinvar"
FILENAME="HROM_yesambi_PGPC-02_recovered_TP_merged.vcf.gz"

REF="/home/zunuan/anaconda3/salfilter/reference/${REFERENCE}/WholeGenomeFasta/genome.fa"
DATA_DIR="/home/zunuan/anaconda3/salfilter/data"
DATA_SOURCES_DIR="/home/zunuan/resource/funcotator_ClinVar_only"

# search for vcf file using *.vcf.idx in directory
IN_FILE=${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${FILENAME}
OUT_FILE=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${FILENAME}

# make output directory if it doesnt exist
if [ ! -d ${OUT_FILE_DIR} ]; then
  mkdir -p ${OUT_FILE_DIR}
fi

# parse parameters

# using the default directory provided by GATK.
# includes the following annotation databases : 
# gnomAD_genome
# gnomAD_exome
# clinvar, gencode, lmm_known, acmg_lof, acmg_rec

gatk Funcotator \
  --variant $IN_FILE \
  --reference $REF \
  --ref-version hg38 \
  --data-sources-path $DATA_SOURCES_DIR \
  --output $OUT_FILE \
  --output-file-format VCF
