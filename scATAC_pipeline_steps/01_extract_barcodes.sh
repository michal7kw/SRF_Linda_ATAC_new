#!/bin/bash
#SBATCH --job-name=extract_barcodes
#SBATCH --output=logs/extract_barcodes_%A_%a.out
#SBATCH --error=logs/extract_barcodes_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"

echo "========================================="
echo "Step 1: Extracting barcodes for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create output directory
mkdir -p "$PROCESSED_DATA_DIR"

# Extract barcodes from RIGHTMOST 16bp of 24bp R2 (positions 9-24)
EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

echo "DEBUG: R2 structure according to facility:"
echo "  [8bp SPACER][16bp 10x barcode] = 24bp total"
echo "  Extracting rightmost 16bp (positions 9-24)"

if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "DEBUG: Extracting rightmost 16bp from R2..."
    
    # Extract rightmost 16bp from 24bp R2 read
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
        awk '{
            if(NR%4==1) print $0;
            else if(NR%4==2) print substr($0,9,16);
            else if(NR%4==3) print $0;
            else if(NR%4==0) print substr($0,9,16);
        }' | gzip > "$EXTRACTED_BC_FILE"
        
    echo "DEBUG: Barcode extraction completed"
else
    echo "DEBUG: Using existing extracted barcodes: $EXTRACTED_BC_FILE"
fi

echo "Extracted barcodes saved to: $EXTRACTED_BC_FILE"

# Quick validation
BC_COUNT=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')
SAMPLE_SEQ=$(zcat "$EXTRACTED_BC_FILE" | awk 'NR==2' | head -1)
SEQ_LENGTH=${#SAMPLE_SEQ}

echo "Validation:"
echo "  Total barcodes: $BC_COUNT"
echo "  Sample barcode: $SAMPLE_SEQ"
echo "  Barcode length: $SEQ_LENGTH bp"

if [[ $SEQ_LENGTH -ne 16 ]]; then
    echo "ERROR: Expected 16bp barcodes, got ${SEQ_LENGTH}bp"
    exit 1
fi

echo "========================================="
echo "Step 1 complete for $SAMPLE"
echo "End time: $(date)"
echo "========================================="