#!/bin/bash
# Trim adapters from ATAC-seq FASTQ files

# Install cutadapt if not present
if ! command -v cutadapt &> /dev/null
then
    echo "cutadapt could not be found, installing..."
    pip install --user cutadapt
fi

# Configuration
DATA_DIR="ATAC_data/nestin"
OUTPUT_DIR="ATAC_data/nestin_trimmed"
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" # Standard Illumina adapter

mkdir -p $OUTPUT_DIR

# Trim adapters for each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Trimming adapters for sample: $SAMPLE"
    ~/.local/bin/cutadapt -a $ADAPTER -o $OUTPUT_DIR/${SAMPLE}_R1_001.fastq.gz $DATA_DIR/${SAMPLE}_R1_001.fastq.gz
    ~/.local/bin/cutadapt -a $ADAPTER -o $OUTPUT_DIR/${SAMPLE}_R2_001.fastq.gz $DATA_DIR/${SAMPLE}_R2_001.fastq.gz
    ~/.local/bin/cutadapt -a $ADAPTER -o $OUTPUT_DIR/${SAMPLE}_R3_001.fastq.gz $DATA_DIR/${SAMPLE}_R3_001.fastq.gz
    cp $DATA_DIR/${SAMPLE}_I1_001.fastq.gz $OUTPUT_DIR/
done

echo "Adapter trimming complete."