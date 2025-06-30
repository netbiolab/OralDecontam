#!/bin/bash
  
REFERENCE=$1   # "1.gatk_GRCh38"
CALLER=$2      # "1.HC"
PREV=$3        # "2-1.funcotator"
NEXT=$4        # "3.bcftools_norm_annot"
SAMPLE=$5      # "raw_PGPC-02"
SOURCE=$6      # "2.buccal"

REF="/home/zunuan/anaconda3/salfilter/reference/${REFERENCE}/WholeGenomeFasta/genome.fa"
DATA_DIR="/home/zunuan/anaconda3/salfilter/data"
SCRIPT_DIR=/home/zunuan/anaconda3/salfilter/script/

# search for vcf file using *.vcf.idx in directory
IN_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${SAMPLE}/${SOURCE}
IN_FILE=$(ls $IN_FILE_DIR | egrep ".vcf.gz$")
ACCESSION="${IN_FILE%.vcf.gz}"

OUT_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}
# make output directory if it doesnt exist
if [ ! -d ${OUT_FILE_DIR} ]; then
  mkdir -p ${OUT_FILE_DIR}
fi

# parse parameters
in_vcf=${IN_FILE_DIR}/${ACCESSION}.vcf.gz
temp_norm_vcf=${OUT_FILE_DIR}/${ACCESSION}_norm.tmp.vcf.gz
bcftools norm \
        -m -both \
        -o ${temp_norm_vcf} \
        ${in_vcf}

out_annot_tab=${OUT_FILE_DIR}/annot.tab
out_annot_hdr=${OUT_FILE_DIR}/annot.hdr
python3 ${SCRIPT_DIR}/0.pipeline/helper_get_annots_tab_from_funcotated_vcf.py \
        ${temp_norm_vcf} \
        ${out_annot_tab} \
	${out_annot_hdr}

bgzip ${out_annot_tab}
tabix -s1 -b2 -e3 \
	${out_annot_tab}.gz

out_annot_vcf=${OUT_FILE_DIR}/${ACCESSION}_annot.vcf.gz

bcftools annotate \
	-a ${out_annot_tab}.gz \
	-h ${out_annot_hdr} \
	-c CHROM,FROM,TO,-,ALT,EXOME_AF_POPMAX,GENOME_AF_POPMAX \
	-O z \
	-o ${out_annot_vcf} \
	${temp_norm_vcf}

rm ${temp_norm_vcf}
