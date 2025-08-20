#!/bin/bash
#SBATCH --job-name=cellranger_arc
#SBATCH --output=logs/cellranger_arc_%A_%a.out
#SBATCH --error=logs/cellranger_arc_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --partition=workq

# CellRanger ARC Processing Script for ARC-v1 Multiome ATAC-only Data
# This script uses a samplesheet to define input files.

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

# Configuration
CELLRANGER_ARC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-arc-2.0.2/bin/cellranger-arc"
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A"
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/nestin_prepared"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_arc_output"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp"
SAMPLESHEET="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/samplesheet.csv"

# Sample names
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing ARC-v1 multiome ATAC-only sample: $SAMPLE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "logs"

# Generate unique timestamp
TIMESTAMP=$(date +%Y%m%d%H%M%S)
UNIQUE_ID="${SAMPLE}_arc_${TIMESTAMP}"
LOCAL_TMP="$TMP_DIR/cellranger_arc_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"

# Create working directories
mkdir -p "$LOCAL_TMP"

# Create libraries CSV for ARC with ATAC-only configuration
echo "Creating libraries CSV for ARC-v1 multiome ATAC-only processing..."
LIBRARIES_CSV="$LOCAL_TMP/${SAMPLE}_libraries.csv"
{
    echo "fastqs,sample,library_type"
    echo "$DATA_DIR,$SAMPLE,Chromatin Accessibility"
} > "$LIBRARIES_CSV"

echo "Libraries CSV content:"
cat "$LIBRARIES_CSV"

# Validate CellRanger ARC executable
if [[ ! -x "$CELLRANGER_ARC" ]]; then
    echo "ERROR: CellRanger ARC not found or not executable: $CELLRANGER_ARC"
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

echo "Starting CellRanger ARC for ARC-v1 multiome ATAC-only data..."
echo "Sample: $SAMPLE"
echo "Reference: $REF"
echo "Libraries CSV: $LIBRARIES_CSV"
echo "Working directory: $LOCAL_TMP"

# Run CellRanger ARC count
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
    
    RESULTS_DIR="$OUTPUT_DIR/${SAMPLE}_arc_results"
    mkdir -p "$RESULTS_DIR"
    
    echo "Copying results to: $RESULTS_DIR"
    if [[ -d "$UNIQUE_ID/outs" ]]; then
        cp -r "$UNIQUE_ID/outs/"* "$RESULTS_DIR/"
        echo "Results successfully copied to $RESULTS_DIR"
    else
        echo "ERROR: Expected output directory not found: $UNIQUE_ID/outs"
        exit 1
    fi
    
else
    echo "ERROR: CellRanger ARC failed for sample: $SAMPLE"
    echo "Check logs for more details"
    if [[ -f "$UNIQUE_ID/_log" ]]; then
        cp "$UNIQUE_ID/_log" "$OUTPUT_DIR/${SAMPLE}_error.log"
    fi
    exit 1
fi

echo "Processing completed for $SAMPLE"