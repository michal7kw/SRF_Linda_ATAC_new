#!/bin/bash
#SBATCH --job-name=trim_adapters_atac
#SBATCH --output=logs/trim_adapters_%A_%a.out
#SBATCH --error=logs/trim_adapters_%A_%a.err
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# Parallelized Adapter Trimming Script for ATAC-seq FASTQ files
# Each FASTQ file is processed as a separate array job for maximum parallelization

set -euo pipefail

# Cleanup function
cleanup() {
    EXIT_CODE=$?
    echo "Performing cleanup..."
    pkill -P $$ || true
    echo "Cleanup complete."
    exit $EXIT_CODE
}
trap cleanup EXIT

# Install cutadapt if not present
if ! command -v cutadapt &> /dev/null; then
    echo "cutadapt could not be found, installing..."
    pip install --user cutadapt
fi

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data"
DATA_DIR="$BASE_DIR/nestin"
OUTPUT_DIR="$BASE_DIR/nestin_trimmed"
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" # Standard Illumina adapter

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Define all files to process (each as separate array task)
# Format: "SAMPLE:READ:ACTION"
FILES=(
    "R26-Nestin-Ctrl-adult:R1:trim"
    "R26-Nestin-Ctrl-adult:R2:trim"  
    "R26-Nestin-Ctrl-adult:R3:trim"
    "R26-Nestin-Ctrl-adult:I1:copy"
    "R26-Nestin-Mut-adult:R1:trim"
    "R26-Nestin-Mut-adult:R2:trim"
    "R26-Nestin-Mut-adult:R3:trim"
    "R26-Nestin-Mut-adult:I1:copy"
)

# Get file info for this array job
FILE_INFO="${FILES[$SLURM_ARRAY_TASK_ID]}"
IFS=':' read -r SAMPLE READ ACTION <<< "$FILE_INFO"

echo "Processing Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample: $SAMPLE"
echo "Read: $READ" 
echo "Action: $ACTION"
echo "Using $SLURM_CPUS_PER_TASK cores"

# Define file paths
INPUT_FILE="$DATA_DIR/${SAMPLE}_${READ}_001.fastq.gz"
OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE}_${READ}_001.fastq.gz"

# Verify input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"

# Process file based on action
if [[ "$ACTION" == "trim" ]]; then
    echo "Trimming adapters from $READ read (genomic DNA)"
    ~/.local/bin/cutadapt \
        --cores="$SLURM_CPUS_PER_TASK" \
        -a "$ADAPTER" \
        -o "$OUTPUT_FILE" \
        "$INPUT_FILE"
    
    if [[ $? -eq 0 ]]; then
        echo "Successfully trimmed: $OUTPUT_FILE"
        echo "Output file size: $(du -h "$OUTPUT_FILE" | cut -f1)"
    else
        echo "ERROR: Failed to trim $INPUT_FILE"
        exit 1
    fi
    
elif [[ "$ACTION" == "copy" ]]; then
    echo "Copying $READ read (cell barcodes) without trimming"
    cp "$INPUT_FILE" "$OUTPUT_FILE"
    
    if [[ $? -eq 0 ]]; then
        echo "Successfully copied: $OUTPUT_FILE"
        echo "Output file size: $(du -h "$OUTPUT_FILE" | cut -f1)"
    else
        echo "ERROR: Failed to copy $INPUT_FILE"
        exit 1
    fi
fi

# Verify output file was created
if [[ -f "$OUTPUT_FILE" ]]; then
    echo "SUCCESS: Processing completed for $SAMPLE $READ"
    echo "Final output: $OUTPUT_FILE"
else
    echo "ERROR: Output file not created: $OUTPUT_FILE"
    exit 1
fi