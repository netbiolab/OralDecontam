#!/bin/bash
### usage : PART OF STANDARD PIPELINE ###

# input directories
DATA_DIR=$1
CURRENT=$2
SAMPLE=$3
SOURCE=$4
NEXT=$5

# constant directories
REF=/home/zunuan/anaconda3/salfilter/reference/1.gatk_GRCh38/WholeGenomeFasta/genome.fa
# resource1
HAPMAP=/home/zunuan/resource/GATK_BUNDLE_CURRENT/hapmap_3.3.hg38.vcf.gz
# resource2
HICONF=/home/zunuan/resource/GATK_BUNDLE_CURRENT/1000G_phase1.snps.high_confidence.hg38.vcf.gz
# resource3
MILLS_INDELS=/home/zunuan/resource/GATK_BUNDLE_CURRENT/Mills_and_1000G_gold_standard.indels.hg38.vcf

# search for bam file in directory
VCFIDX_DIR=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}
VCFIDX=$(ls $VCFIDX_DIR | egrep ".vcf.gz$")
ACCESSION="${VCFIDX%.vcf.gz}"

# make output directory if it doesnt exist
#if [ ! -d ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE} ]; then
#  mkdir -p ${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}
#fi

# parse parameters
in_vcf=${DATA_DIR}/${CURRENT}/${SAMPLE}/${SOURCE}/${ACCESSION}.vcf.gz
inter_vcf=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}.tmp.vcf
inter_vcf_idx=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}.tmp.vcf.idx
out_filtered_vcf=${DATA_DIR}/${NEXT}/${SAMPLE}/${SOURCE}/${ACCESSION}_nofilt.vcf


# 1D model with pre-trained architecture
gatk CNNScoreVariants \
	-V $in_vcf \
	-R $REF \
	-O $inter_vcf

gatk FilterVariantTranches \
	-V $inter_vcf \
	--resource $MILLS_INDELS \
	--resource $HICONF \
	--resource $HAPMAP \
	--info-key CNN_1D \
  --snp-tranche 99.95 \
	--indel-tranche 99.4 \
	-O $out_filtered_vcf
# --invalidate-previous-filters : remove existing filters from the input VCF.

bgzip $out_filtered_vcf
tabix -p vcf ${DATA_DIR}.gz 

# remove intermediate files
rm $inter_vcf
rm $inter_vcf_idx
