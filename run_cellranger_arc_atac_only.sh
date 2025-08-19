#!/bin/bash
#SBATCH --job-name=cellranger_arc_atac_only
#SBATCH --output=logs/cellranger_arc_atac_only_%A_%a.out
#SBATCH --error=logs/cellranger_arc_atac_only_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --partition=workq

# CellRanger ARC Processing Script for ARC-v1 Multiome ATAC-only Data
# Based on successful RNA processing pattern that uses --chemistry=ARC-v1

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

# Environment variables for CellRanger (from RNA script)
export CELLRANGER_COPY_MODE=copy
export CELLRANGER_USE_HARDLINKS=false

# Configuration - Use CellRanger ARC for ARC-v1 multiome data
CELLRANGER_ARC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-arc-2.0.2/bin/cellranger-arc"
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A"
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/nestin"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_arc_atac_only_output"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp"

# Sample names (same as RNA processing)
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing ARC-v1 multiome ATAC-only sample: $SAMPLE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "logs"

# Generate unique timestamp
TIMESTAMP=$(date +%Y%m%d%H%M%S)
UNIQUE_ID="${SAMPLE}_arc_atac_only_${TIMESTAMP}"
FASTQ_DIR="$OUTPUT_DIR/${SAMPLE}_fastq"
LOCAL_TMP="$TMP_DIR/cellranger_arc_atac_only_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"
LIBRARIES_CSV="$LOCAL_TMP/${SAMPLE}_libraries.csv"

# Create working directories
mkdir -p "$FASTQ_DIR" "$LOCAL_TMP"

echo "Creating symbolic links for ARC-v1 multiome ATAC files..."

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

# Create symbolic links with proper 10x naming convention
# Following the pattern from RNA script but for ATAC reads
for READ in "I1" "R1" "R2" "R3"; do
    SOURCE_FILE="$DATA_DIR/${SAMPLE}_${READ}_001.fastq.gz"
    TARGET_FILE="$FASTQ_DIR/${SAMPLE}_S1_L001_${READ}_001.fastq.gz"
    ln -sf "$SOURCE_FILE" "$TARGET_FILE"
done

# Create libraries CSV for ARC with ATAC-only configuration
# This is the key difference - we specify only Chromatin Accessibility
echo "Creating libraries CSV for ARC-v1 multiome ATAC-only processing..."
cat > "$LIBRARIES_CSV" <<EOF
fastqs,sample,library_type
$FASTQ_DIR,$SAMPLE,Chromatin Accessibility
EOF

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

# Calculate memory (matching RNA script pattern)
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
# Note: We DON'T use --chemistry flag here as it's not supported by cellranger-arc
# The ARC-v1 chemistry is automatically detected from the read structure
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
    
    # Create results directory
    RESULTS_DIR="$OUTPUT_DIR/${SAMPLE}_arc_atac_only_results"
    mkdir -p "$RESULTS_DIR"
    
    # Copy results
    echo "Copying results to: $RESULTS_DIR"
    if [[ -d "$UNIQUE_ID/outs" ]]; then
        # Copy all essential outputs
        cp -r "$UNIQUE_ID/outs/"* "$RESULTS_DIR/"
        echo "Results successfully copied to $RESULTS_DIR"
        
        # List key ATAC output files
        echo "Key ATAC output files created:"
        for file in "web_summary.html" "summary.csv" "atac_fragments.tsv.gz" "atac_peaks.bed" "atac_peak_annotation.tsv"; do
            if [[ -f "$RESULTS_DIR/$file" ]]; then
                echo "  ✓ $file"
            else
                echo "  ✗ $file (may not be present for ATAC-only)"
            fi
        done
        
        # Also check for the filtered feature matrix
        if [[ -d "$RESULTS_DIR/filtered_feature_bc_matrix" ]]; then
            echo "  ✓ filtered_feature_bc_matrix/"
        elif [[ -f "$RESULTS_DIR/filtered_feature_bc_matrix.h5" ]]; then
            echo "  ✓ filtered_feature_bc_matrix.h5"
        fi
    else
        echo "ERROR: Expected output directory not found: $UNIQUE_ID/outs"
        exit 1
    fi
    
else
    echo "ERROR: CellRanger ARC failed for sample: $SAMPLE"
    echo "Check logs for more details"
    # Save the error log
    if [[ -f "$UNIQUE_ID/_log" ]]; then
        cp "$UNIQUE_ID/_log" "$OUTPUT_DIR/${SAMPLE}_error.log"
    fi
    exit 1
fi

echo "Processing completed for $SAMPLE"