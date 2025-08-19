#!/bin/bash
#SBATCH --job-name=cellranger_arc_mock_rna
#SBATCH --output=logs/cellranger_arc_mock_rna_%A_%a.out
#SBATCH --error=logs/cellranger_arc_mock_rna_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --partition=workq

# CellRanger ARC Processing with Mock RNA Libraries
# This script creates minimal mock RNA files to satisfy CellRanger ARC requirements
# for processing ARC-v1 multiome ATAC-only data

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

# Environment variables
export CELLRANGER_COPY_MODE=copy
export CELLRANGER_USE_HARDLINKS=false

# Configuration
CELLRANGER_ARC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-arc-2.0.2/bin/cellranger-arc"
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A"
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/nestin"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_arc_mock_output"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp"

# Sample names
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing ARC-v1 multiome sample with mock RNA: $SAMPLE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "logs"

# Generate unique timestamp
TIMESTAMP=$(date +%Y%m%d%H%M%S)
UNIQUE_ID="${SAMPLE}_arc_mock_${TIMESTAMP}"
ATAC_FASTQ_DIR="$OUTPUT_DIR/${SAMPLE}_atac_fastq"
RNA_FASTQ_DIR="$OUTPUT_DIR/${SAMPLE}_rna_fastq"
LOCAL_TMP="$TMP_DIR/cellranger_arc_mock_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"
LIBRARIES_CSV="$LOCAL_TMP/${SAMPLE}_libraries.csv"

# Create working directories
mkdir -p "$ATAC_FASTQ_DIR" "$RNA_FASTQ_DIR" "$LOCAL_TMP"

echo "Setting up ATAC files..."
# Verify and link ATAC files
for READ in "I1" "R1" "R2" "R3"; do
    SOURCE_FILE="$DATA_DIR/${SAMPLE}_${READ}_001.fastq.gz"
    if [[ -f "$SOURCE_FILE" ]]; then
        TARGET_FILE="$ATAC_FASTQ_DIR/${SAMPLE}_S1_L001_${READ}_001.fastq.gz"
        ln -sf "$SOURCE_FILE" "$TARGET_FILE"
        echo "  Linked: $READ"
    else
        echo "  ERROR: Missing $SOURCE_FILE"
        exit 1
    fi
done

echo "Creating minimal mock RNA files..."
# Create minimal valid FASTQ files for RNA
# These will have the same number of reads as ATAC but with dummy sequences
# This allows CellRanger ARC to run and process the ATAC data

# Function to create a minimal gzipped FASTQ with dummy data
create_mock_fastq() {
    local output_file=$1
    local read_type=$2
    
    # Create a minimal FASTQ with a few dummy reads
    (
        # Add 100 dummy reads to satisfy minimum requirements
        for i in {1..100}; do
            echo "@DUMMY:1:DUMMY:1:1:1:1 1:N:0:DUMMY"
            if [[ "$read_type" == "R1" ]]; then
                # R1 for RNA should be 28bp in ARC-v1
                echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            elif [[ "$read_type" == "R2" ]]; then
                # R2 for RNA should be 90bp in ARC-v1
                echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            elif [[ "$read_type" == "I1" ]]; then
                # I1 is 8bp
                echo "NNNNNNNN"
            fi
            echo "+"
            if [[ "$read_type" == "R1" ]]; then
                echo "IIIIIIIIIIIIIIIIIIIIIIIIIIII"
            elif [[ "$read_type" == "R2" ]]; then
                echo "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
            elif [[ "$read_type" == "I1" ]]; then
                echo "IIIIIIII"
            fi
        done
    ) | gzip > "$output_file"
}

# Create mock RNA FASTQ files
for READ in "I1" "R1" "R2"; do
    MOCK_FILE="$RNA_FASTQ_DIR/${SAMPLE}_S1_L001_${READ}_001.fastq.gz"
    create_mock_fastq "$MOCK_FILE" "$READ"
    echo "  Created mock RNA $READ"
done

# Create libraries CSV with both ATAC and mock RNA
echo "Creating libraries CSV..."
cat > "$LIBRARIES_CSV" <<EOF
fastqs,sample,library_type
$ATAC_FASTQ_DIR,$SAMPLE,Chromatin Accessibility
$RNA_FASTQ_DIR,$SAMPLE,Gene Expression
EOF

echo "Libraries CSV content:"
cat "$LIBRARIES_CSV"

# Validate executables and reference
if [[ ! -x "$CELLRANGER_ARC" ]]; then
    echo "ERROR: CellRanger ARC not found: $CELLRANGER_ARC"
    exit 1
fi

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

echo "Starting CellRanger ARC with mock RNA libraries..."
echo "Sample: $SAMPLE"
echo "Reference: $REF"
echo "Libraries CSV: $LIBRARIES_CSV"

# Run CellRanger ARC count
# The mock RNA data will fail QC but ATAC should process correctly
"$CELLRANGER_ARC" count \
    --id="$UNIQUE_ID" \
    --reference="$REF" \
    --libraries="$LIBRARIES_CSV" \
    --localcores="$SLURM_CPUS_PER_TASK" \
    --localmem="$CELLRANGER_MEM_GB" \
    --min-atac-count=100 \
    --min-gex-count=10 \
    --disable-ui

# Check if processing succeeded
if [[ $? -eq 0 ]]; then
    echo "CellRanger ARC processing completed for sample: $SAMPLE"
    
    # Create results directory
    RESULTS_DIR="$OUTPUT_DIR/${SAMPLE}_results"
    mkdir -p "$RESULTS_DIR"
    
    # Copy ATAC results only (RNA will be garbage)
    echo "Copying ATAC results to: $RESULTS_DIR"
    if [[ -d "$UNIQUE_ID/outs" ]]; then
        # Copy ATAC-specific outputs
        for file in "atac_fragments.tsv.gz" "atac_fragments.tsv.gz.tbi" \
                   "atac_peaks.bed" "atac_peak_annotation.tsv" \
                   "web_summary.html" "summary.csv"; do
            if [[ -f "$UNIQUE_ID/outs/$file" ]]; then
                cp "$UNIQUE_ID/outs/$file" "$RESULTS_DIR/"
                echo "  ✓ Copied $file"
            fi
        done
        
        # Copy filtered matrix if exists
        if [[ -d "$UNIQUE_ID/outs/filtered_feature_bc_matrix" ]]; then
            cp -r "$UNIQUE_ID/outs/filtered_feature_bc_matrix" "$RESULTS_DIR/"
            echo "  ✓ Copied filtered_feature_bc_matrix/"
        fi
    else
        echo "ERROR: Output directory not found: $UNIQUE_ID/outs"
        exit 1
    fi
else
    echo "ERROR: CellRanger ARC failed for sample: $SAMPLE"
    exit 1
fi

echo "Processing completed for $SAMPLE (ATAC data extracted)"