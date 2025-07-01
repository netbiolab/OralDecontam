#!/bin/bash

REFERENCE="2.no_alt_GRCh38"
CALLER="7.DV_uncommon"
PREV="3.bcftools_norm_annot"
NEXT=""

REF="/home/zunuan/anaconda3/salfilter/reference/${REFERENCE}/WholeGenomeFasta/genome.fa"
DATA_DIR="/home/zunuan/anaconda3/salfilter/data"

# First set of samples and sources
SAMPLE_LIST1=(
    eHOMD_yesambi_PGPC-02
    eHOMD_yesambi_PGPC-05
    eHOMD_yesambi_PGPC-06
    eHOMD_yesambi_PGPC-50
    HROM_yesambi_PGPC-02
    HROM_yesambi_PGPC-05
    HROM_yesambi_PGPC-06
    HROM_yesambi_PGPC-50
)
SOURCE_LIST1=(
    1.saliva
    2.buccal
    3.saliva_enriched
    4.buccal_enriched
)

# Second set of samples and sources
SAMPLE_LIST2=(
    raw_PGPC-02
    raw_PGPC-05
    raw_PGPC-06
    raw_PGPC-50
)
SOURCE_LIST2=(
    0.blood
    1.saliva
    2.buccal
    3.saliva_enriched
    4.buccal_enriched
)

# Function to process each SAMPLE-SOURCE pair
process_sample() {
    SAMPLE=$1
    SOURCE=$2

    echo "Processing SAMPLE: $SAMPLE | SOURCE: $SOURCE"

    # Input
    IN_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${PREV}/${SAMPLE}/${SOURCE}
    IN_FILE=$(ls $IN_FILE_DIR | egrep ".vcf.gz$")
    ACCESSION="${IN_FILE%.vcf.gz}"
    in_vcf=${IN_FILE_DIR}/${ACCESSION}.vcf.gz

    # Common variants (MAF ≥ 0.05)
    echo "Filtering common variants (MAF ≥ 0.05)"
    NEXT="3-1.common"
    OUT_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}
    mkdir -p ${OUT_FILE_DIR}
    out_vcf=${OUT_FILE_DIR}/${ACCESSION}.vcf.gz

    bcftools filter \
        -i "INFO/EXOME_AF_POPMAX>=0.05 || INFO/GENOME_AF_POPMAX>=0.05" \
        ${in_vcf} -Oz -o ${out_vcf}
    tabix -p vcf ${out_vcf}
    COMMON_COUNT=$(bcftools view -H ${out_vcf} | wc -l)

    # Uncommon variants (0.01 ≤ MAF < 0.05)
    echo "Filtering uncommon variants (0.01 ≤ MAF < 0.05)"
    NEXT="3-2.uncommon"
    OUT_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}
    mkdir -p ${OUT_FILE_DIR}
    temp_vcf=${OUT_FILE_DIR}/${ACCESSION}.tmp.vcf.gz
    out_vcf=${OUT_FILE_DIR}/${ACCESSION}.vcf.gz

    bcftools filter \
        -e "INFO/EXOME_AF_POPMAX>=0.05 || INFO/GENOME_AF_POPMAX>=0.05" \
        ${in_vcf} -Oz -o ${temp_vcf}
    bcftools filter \
        -e "INFO/EXOME_AF_POPMAX<0.01 && INFO/GENOME_AF_POPMAX<0.01" \
        ${temp_vcf} -Oz -o ${out_vcf}
    tabix -p vcf ${out_vcf}
    rm ${temp_vcf}
    UNCOMMON_COUNT=$(bcftools view -H ${out_vcf} | wc -l)

    # Rare variants (MAF < 0.01)
    echo "Filtering rare variants (MAF < 0.01)"
    NEXT="3-3.rare"
    OUT_FILE_DIR=${DATA_DIR}/${REFERENCE}/${CALLER}/${NEXT}/${SAMPLE}/${SOURCE}
    mkdir -p ${OUT_FILE_DIR}
    out_vcf=${OUT_FILE_DIR}/${ACCESSION}.vcf.gz

    bcftools filter \
        -i "INFO/EXOME_AF_POPMAX<0.01 && INFO/GENOME_AF_POPMAX<0.01" \
        ${in_vcf} -Oz -o ${out_vcf}
    tabix -p vcf ${out_vcf}
    RARE_COUNT=$(bcftools view -H ${out_vcf} | wc -l)

    # Sanity check
    TOTAL_COUNT=$(bcftools view -H ${in_vcf} | wc -l)
    if [[ $((COMMON_COUNT + UNCOMMON_COUNT + RARE_COUNT)) -eq $TOTAL_COUNT ]]; then
        echo "Variant count sanity check SUCCESS"
    else
        echo "ERROR: The sum of common, uncommon, and rare variants does not match the total number of variants."
    fi
}

# Loop through first set of samples and sources
for SAMPLE in "${SAMPLE_LIST1[@]}"; do
    for SOURCE in "${SOURCE_LIST1[@]}"; do
        process_sample "$SAMPLE" "$SOURCE"
    done
done

# Loop through second set of samples and sources
for SAMPLE in "${SAMPLE_LIST2[@]}"; do
    for SOURCE in "${SOURCE_LIST2[@]}"; do
        process_sample "$SAMPLE" "$SOURCE"
    done
done
