#!/bin/bash
#SBATCH --job-name=cellranger_atac_arcv1
#SBATCH --output=logs/cellranger_atac_arcv1_%A_%a.out
#SBATCH --error=logs/cellranger_atac_arcv1_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --partition=workq

# CellRanger ATAC Processing Script for ARC-v1 Chemistry
# Based on successful RNA processing pattern from SRF_Linda_RNA

set -euo pipefail

# Cleanup function
cleanup() {
    EXIT_CODE=$?
    echo "Performing cleanup..."
    pkill -P $$ || true
    if [[ $EXIT_CODE -eq 0 ]]; then
        echo "Job finished successfully. Removing temporary directory: $LOCAL_TMP"
        rm -rf "$LOCAL_TMP"
    else
        echo "Job failed with exit code $EXIT_CODE. Temporary directory preserved for inspection: $LOCAL_TMP"
    fi
    echo "Cleanup complete."
}
trap cleanup EXIT

# Environment variables for CellRanger
export CELLRANGER_COPY_MODE=copy
export CELLRANGER_USE_HARDLINKS=false

# Configuration - following successful RNA processing pattern
CELLRANGER_ARC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-arc-2.0.2/cellranger-arc"
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A"
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_atac_arcv1_output_fixed"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp_fixed"

# Sample names (matching RNA processing pattern)
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
LOCAL_TMP="$TMP_DIR/cellranger_atac_arcv1_${SAMPLE}_${SLURM_ARRAY_TASK_ID}" # Define LOCAL_TMP earlier

# Barcode whitelist for ARC-v1 chemistry
BARCODE_WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"
BARCODE_WHITELIST_FILE="atac_737K-arc-v1.txt"
BARCODE_WHITELIST="$BARCODE_WHITELIST_DIR/$BARCODE_WHITELIST_FILE"

mkdir -p "$BARCODE_WHITELIST_DIR"

if [[ ! -f "$BARCODE_WHITELIST" ]]; then
    echo "Downloading ARC-v1 ATAC barcode whitelist..."
    wget -O "$BARCODE_WHITELIST_DIR/${BARCODE_WHITELIST_FILE}.gz" https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz
    gunzip "$BARCODE_WHITELIST_DIR/${BARCODE_WHITELIST_FILE}.gz"
    echo "Barcode whitelist downloaded and unzipped to $BARCODE_WHITELIST"
else
    echo "Barcode whitelist already exists: $BARCODE_WHITELIST"
fi

# Sample names are now defined earlier
echo "Processing ATAC sample with ARC-v1 chemistry: $SAMPLE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

# Generate unique timestamp
TIMESTAMP=$(date +%Y%m%d%H%M%S)
UNIQUE_ID="${SAMPLE}_atac_arcv1_${TIMESTAMP}"
FASTQ_DIR="$OUTPUT_DIR/${SAMPLE}_fastq"
# LOCAL_TMP is now defined earlier

# Create working directories
mkdir -p "$FASTQ_DIR" "$LOCAL_TMP"

echo "Creating symbolic links for ARC-v1 ATAC files..."

# Verify input files exist
echo "Verifying input files exist:"
for READ in "I1" "R1" "R2" "R3"; do
    SOURCE_FILE="$DATA_DIR/${SAMPLE}_${READ}_001.fastq.gz"
    if [[ -f "$SOURCE_FILE" ]]; then
        echo "  Found: $SOURCE_FILE"
    else
        echo "  ERROR: Missing $SOURCE_FILE"
        exit 1
    fi
done

# Create symbolic links following RNA processing pattern
# The key insight: this is ARC-v1 multiome data, so file structure should match
# Create symbolic links following RNA processing pattern
# The key insight: this is ARC-v1 multiome data, so file structure should match
for READ in "I1" "R1" "R3"; do # R2 is handled separately
    SOURCE_FILE="$DATA_DIR/${SAMPLE}_${READ}_001.fastq.gz"
    TARGET_FILE="$FASTQ_DIR/${SAMPLE}_S1_L001_${READ}_001.fastq.gz"
    ln -sf "$SOURCE_FILE" "$TARGET_FILE"
done

# For R2, use the 16bp extracted barcode file
SOURCE_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_16bp.fastq.gz"
TARGET_FILE="$FASTQ_DIR/${SAMPLE}_S1_L001_R2_001.fastq.gz"
ln -sf "$SOURCE_FILE" "$TARGET_FILE"

# Validate CellRanger ATAC executable
if [[ ! -x "$CELLRANGER_ARC" ]]; then
    echo "ERROR: CellRanger ARC not found: $CELLRANGER_ARC"
    exit 1
fi

# Validate reference genome
if [[ ! -d "$REF" ]]; then
    echo "ERROR: Reference genome not found: $REF"
    exit 1
fi

# Calculate memory
TOTAL_MEM_GB=$(($SLURM_MEM_PER_NODE / 1024))
CELLRANGER_MEM_GB=$((TOTAL_MEM_GB - 16))

echo "Using $SLURM_CPUS_PER_TASK cores and ${CELLRANGER_MEM_GB}GB memory"

# Change to temporary directory
cd "$LOCAL_TMP"

# Clean up any existing pipestance
if [[ -d "$UNIQUE_ID" ]]; then
    echo "Removing existing pipestance directory: $UNIQUE_ID"
    rm -rf "$UNIQUE_ID"
fi

echo "Starting CellRanger ATAC with ARC-v1 chemistry compatibility..."
echo "Sample: $SAMPLE"
echo "Reference: $REF"
echo "FASTQ directory: $FASTQ_DIR"
echo "Working directory: $LOCAL_TMP"

# Try CellRanger ATAC with explicit chemistry and parameters optimized for ARC-v1
# Following the pattern from successful RNA processing
LIBRARIES_CSV="$LOCAL_TMP/libraries.csv"
echo "fastqs,sample,library_type" > "$LIBRARIES_CSV"
echo "$FASTQ_DIR,$SAMPLE,\"Chromatin Accessibility\"" >> "$LIBRARIES_CSV"
echo "$FASTQ_DIR,$SAMPLE,\"Gene Expression\"" >> "$LIBRARIES_CSV"

echo "Generated libraries.csv:"
cat "$LIBRARIES_CSV"

"$CELLRANGER_ARC" count \
    --id="$UNIQUE_ID" \
    --reference="$REF" \
    --libraries="$LIBRARIES_CSV" \
    --localcores="$SLURM_CPUS_PER_TASK" \
    --localmem="$CELLRANGER_MEM_GB" \
    --disable-ui

# Check if processing succeeded
if [[ $? -eq 0 ]]; then
    echo "CellRanger ARC processing completed successfully for sample: $SAMPLE"
    
    # Create results directory following RNA pattern
    ATAC_DIR="$OUTPUT_DIR/${SAMPLE}_atac_arcv1_results"
    mkdir -p "$ATAC_DIR"
    
    # Copy results
    echo "Copying results to: $ATAC_DIR"
    if [[ -d "$UNIQUE_ID/outs" ]]; then
        # Copy essential ATAC outputs
        cp -r "$UNIQUE_ID/outs/"* "$ATAC_DIR/"
        echo "Results successfully copied to $ATAC_DIR"
        
        # Verify key ATAC output files
        echo "Key ATAC output files created:"
        for file in "web_summary.html" "summary.csv" "filtered_peak_bc_matrix.h5" "fragments.tsv.gz" "peaks.bed" "singlecell.csv"; do
            if [[ -f "$ATAC_DIR/$file" ]]; then
                echo "  ✓ $file"
            else
                echo "  ✗ $file (missing)"
            fi
        done
        
        # Copy full output for reference
        rsync -av "$UNIQUE_ID/" "$OUTPUT_DIR/$UNIQUE_ID/"
    else
        echo "ERROR: Expected output directory not found: $UNIQUE_ID/outs"
        exit 1
    fi
    
else
    echo "ERROR: CellRanger ARC failed for sample: $SAMPLE"
    echo "Check logs for more details"
    exit 1
fi

echo "ATAC processing completed for $SAMPLE using ARC-v1 compatible approach"