# OralDecontam

This repository provides an end-to-end pipeline that classifies paired-end FASTQ reads using [Kraken2](https://ccb.jhu.edu/software/kraken2/) and separates the classified reads into three categories:
- **Human (GRCh38)**
- **Non-Human**
- **Ambiguous/Root** (reads only classified to the root taxonomic level)

The tool is built in Python and designed to work with large-scale metagenomic or sequencing data workflows.

## Use Case

Ideal for oral samples before alignment for genotpying purposes.

---

## Requirements

- Python 3.7+
- [`Kraken2`](https://github.com/DerrickWood/kraken2) installed and in `$PATH`
- Optional: `wget`, `curl`, or `aria2c` for DB download

---

## Installation

Clone the repository:

```bash
git clone https://github.com/your-username/kraken2-human-nonhuman-split.git
cd kraken2-human-nonhuman-split
Install any missing Python dependencies:

pip install -r requirements.txt  # if you modularize requirements
Usage
Basic CLI Example
python kraken2_split_pipeline.py \
    --forward sample_R1.fastq.gz \
    --reverse sample_R2.fastq.gz \
    --db /path/to/kraken2_db \
    --prefix Sample01
With Custom Output Prefixes and Directory
python kraken2_split_pipeline.py \
    --forward sample_R1.fastq.gz \
    --reverse sample_R2.fastq.gz \
    --db ./kraken_db \
    --classified classified_output \
    --unclassified unclassified_output \
    --prefix ProjectX \
    --read_dir /data/reads
Output Files
All output files will be written to the current working directory (or an optional output directory, if added):

Kraken2 Outputs
kraken2_output.txt: Kraken2 classification log
kraken2_report.txt: Classification summary
classified_1.fastq, classified_2.fastq: Reads Kraken2 could classify
unclassified_1.fastq, unclassified_2.fastq: Reads Kraken2 couldn't classify
Separated FASTQ Files
Category	Forward Read File	Reverse Read File
Human	Sample01.grch38_1.fastq	Sample01.grch38_2.fastq
Non-Human	Sample01.non-human_1.fastq	Sample01.non-human_2.fastq
Ambiguous/Root	Sample01.ambiguous_1.fastq	Sample01.ambiguous_2.fastq

Taxonomic Classification Criteria
Human reads are identified based on the following markers:

d__Eukaryota, p__Chordata, c__Mammalia, o__Primates,
f__Hominidae, g__Homo, s__Homosapiens_GrCH38.fna, GRCH38.fna
Reads classified only to root (taxid 1) are considered ambiguous.

All remaining reads are treated as non-human.

📚 References
Kraken2: https://ccb.jhu.edu/software/kraken2/
Human Genome GRCh38: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/
