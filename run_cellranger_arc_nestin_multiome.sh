#!/bin/bash
#SBATCH --job-name=cellranger_arc_nestin
#SBATCH --output=logs/cellranger_arc_%A_%a.out
#SBATCH --error=logs/cellranger_arc_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --partition=workq

# CellRanger ARC Processing Script for Nestin Multiome ATAC-only Data
# This script processes multiome data with ATAC-only libraries using CellRanger ARC

set -euo pipefail

# Configuration
CELLRANGER_ARC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-arc-2.0.2/bin/cellranger-arc"
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A"
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/nestin"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_arc_output"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp"
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing multiome ATAC sample: $SAMPLE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

# Generate unique timestamp for this run
TIMESTAMP=$(date +%Y%m%d%H%M%S)
UNIQUE_ID="${SAMPLE}_multiome_${TIMESTAMP}"
FASTQ_DIR="$OUTPUT_DIR/${SAMPLE}_fastq"
LOCAL_TMP="$TMP_DIR/cellranger_arc_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"
LIBRARIES_CSV="$LOCAL_TMP/${SAMPLE}_libraries.csv"

# Cleanup function
cleanup() {
    echo "Cleaning up processes..."
    pkill -P $$ || true
    if [[ -d "$LOCAL_TMP" ]]; then
        rm -rf "$LOCAL_TMP"
    fi
    echo "Cleanup complete."
}
trap cleanup EXIT

# Create temporary working directory
mkdir -p "$LOCAL_TMP"

# Create FASTQ directory and symbolic links
echo "Creating symbolic links for multiome ATAC files..."
mkdir -p "$FASTQ_DIR"

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
for READ in "I1" "R1" "R2" "R3"; do
    SOURCE_FILE="$DATA_DIR/${SAMPLE}_${READ}_001.fastq.gz"
    TARGET_FILE="$FASTQ_DIR/${SAMPLE}_S1_L001_${READ}_001.fastq.gz"
    ln -sf "$SOURCE_FILE" "$TARGET_FILE"
done

# Create libraries CSV for multiome ATAC-only processing
echo "Creating libraries CSV for multiome ATAC processing..."
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

# Calculate memory (leave some headroom)
TOTAL_MEM_GB=$(($SLURM_MEM_PER_NODE / 1024))
CELLRANGER_MEM_GB=$((TOTAL_MEM_GB - 16))

echo "Using $SLURM_CPUS_PER_TASK cores and ${CELLRANGER_MEM_GB}GB memory for Cell Ranger ARC"

# Change to temporary directory for processing
cd "$LOCAL_TMP"

echo "Starting Cell Ranger ARC multiome processing..."
echo "Sample: $SAMPLE"
echo "Reference: $REF"
echo "Libraries CSV: $LIBRARIES_CSV"
echo "Output directory: $LOCAL_TMP"

# Run CellRanger ARC with multiome ATAC-only libraries
"$CELLRANGER_ARC" count \
    --id="$UNIQUE_ID" \
    --reference="$REF" \
    --libraries="$LIBRARIES_CSV" \
    --localcores="$SLURM_CPUS_PER_TASK" \
    --localmem="$CELLRANGER_MEM_GB" \
    --disable-ui

# Check if processing succeeded
if [[ $? -eq 0 ]]; then
    echo "Cell Ranger ARC processing completed successfully for sample: $SAMPLE"
    
    # Create results directory
    RESULTS_DIR="$OUTPUT_DIR/${SAMPLE}_multiome_results"
    mkdir -p "$RESULTS_DIR"
    
    # Copy key output files
    echo "Copying results to: $RESULTS_DIR"
    if [[ -d "$UNIQUE_ID/outs" ]]; then
        # Copy all essential outputs
        cp -r "$UNIQUE_ID/outs/"* "$RESULTS_DIR/"
        echo "Results successfully copied to $RESULTS_DIR"
        
        # List key output files
        echo "Key output files created:"
        for file in "web_summary.html" "summary.csv" "filtered_feature_bc_matrix.h5" "atac_fragments.tsv.gz" "atac_peaks.bed"; do
            if [[ -f "$RESULTS_DIR/$file" ]]; then
                echo "  ✓ $file"
            else
                echo "  ✗ $file (missing)"
            fi
        done
    else
        echo "ERROR: Expected output directory not found: $UNIQUE_ID/outs"
        exit 1
    fi
    
else
    echo "ERROR: Cell Ranger ARC failed for sample: $SAMPLE"
    echo "Check logs for more details"
    exit 1
fi

echo "Processing completed for $SAMPLE"