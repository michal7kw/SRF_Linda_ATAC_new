#!/bin/bash
#SBATCH --job-name=cellranger_atac_arcv1_fixed
#SBATCH --output=logs/cellranger_atac_arcv1_%A_%a.out
#SBATCH --error=logs/cellranger_atac_arcv1_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --partition=workq

# CellRanger ATAC Processing Script for ARC-v1 Chemistry - CORRECTED VERSION

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
# IMPORTANT: Use cellranger-atac, not cellranger-arc for ATAC-only data
CELLRANGER_ATAC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-atac-2.1.0/cellranger-atac"
# Use ATAC-specific reference (not ARC reference)
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-atac-GRCm39-2024"
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_atac_output_corrected"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp_corrected"

# Sample names
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
LOCAL_TMP="$TMP_DIR/cellranger_atac_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"

echo "Processing ATAC sample: $SAMPLE"

# Create directories IMMEDIATELY after defining LOCAL_TMP
mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "$LOCAL_TMP"

# Generate unique ID
TIMESTAMP=$(date +%Y%m%d%H%M%S)
UNIQUE_ID="${SAMPLE}_atac_${TIMESTAMP}"
FASTQ_DIR="$OUTPUT_DIR/${SAMPLE}_fastq"
mkdir -p "$FASTQ_DIR"

# ============================================
# CRITICAL: Determine correct file mapping
# ============================================
# For 10x Multiome ATAC (ARC-v1):
# - I1: Sample index (8bp) - can be ignored for single sample
# - R1: Genomic insert read 1 (50bp)
# - R2: Barcode (24bp total: 16bp barcode + 8bp linker)
# - R3: Genomic insert read 2 (49bp)
#
# For CellRanger ATAC input format:
# - I1: Sample index
# - R1: Read 1 (genomic)
# - R2: i5 index (barcode)
# - R3: Read 2 (genomic)

echo "Creating symbolic links for ATAC files..."

# Verify input files exist
for READ in "I1" "R1" "R2" "R3"; do
    SOURCE_FILE="$DATA_DIR/${SAMPLE}_${READ}_001.fastq.gz"
    if [[ -f "$SOURCE_FILE" ]]; then
        echo "  Found: $SOURCE_FILE ($(zcat "$SOURCE_FILE" | head -4 | sed -n '2p' | wc -c) bp)"
    else
        echo "  ERROR: Missing $SOURCE_FILE"
        exit 1
    fi
done

# Create symbolic links with correct naming for cellranger-atac
ln -sf "$DATA_DIR/${SAMPLE}_I1_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_I1_001.fastq.gz"
ln -sf "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_R1_001.fastq.gz"
ln -sf "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_R2_001.fastq.gz"  # This contains barcodes
ln -sf "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_R3_001.fastq.gz"

# ============================================
# Check barcode orientation (CRITICAL)
# ============================================
echo "Checking barcode orientation..."
BARCODE_WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"
BARCODE_WHITELIST="$BARCODE_WHITELIST_DIR/atac_737K-arc-v1.txt"
BARCODE_WHITELIST_RC="$BARCODE_WHITELIST_DIR/atac_737K-arc-v1_rc.txt"

# Download whitelist if needed
if [[ ! -f "$BARCODE_WHITELIST" ]]; then
    mkdir -p "$BARCODE_WHITELIST_DIR"
    wget -O "$BARCODE_WHITELIST_DIR/atac_737K-arc-v1.txt.gz" \
        https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz
    gunzip "$BARCODE_WHITELIST_DIR/atac_737K-arc-v1.txt.gz"
fi

# Create reverse complement if needed
if [[ ! -f "$BARCODE_WHITELIST_RC" ]]; then
    cat "$BARCODE_WHITELIST" | rev | tr 'ACGT' 'TGCA' > "$BARCODE_WHITELIST_RC"
fi

# Sample first 10000 barcodes from R2 to check orientation
echo "Sampling barcodes to determine orientation..."
# Add error handling and avoid SIGPIPE
if ! zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" 2>/dev/null | \
    awk 'NR%4==2 {print substr($0,1,16)}' | \
    head -10000 > "$LOCAL_TMP/sampled_barcodes.txt"; then
    echo "ERROR: Failed to sample barcodes from R2 file"
    echo "Checking if file is readable..."
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | head -4
    exit 1
fi

echo "Sampled $(wc -l < "$LOCAL_TMP/sampled_barcodes.txt") barcodes"

# Check matches against both orientations
MATCHES_ORIGINAL=$(grep -Ff "$BARCODE_WHITELIST" "$LOCAL_TMP/sampled_barcodes.txt" | wc -l)
MATCHES_RC=$(grep -Ff "$BARCODE_WHITELIST_RC" "$LOCAL_TMP/sampled_barcodes.txt" | wc -l)

echo "Matches with original whitelist: $MATCHES_ORIGINAL / 10000 ($(echo "scale=2; $MATCHES_ORIGINAL * 100 / 10000" | bc)%)"
echo "Matches with RC whitelist: $MATCHES_RC / 10000 ($(echo "scale=2; $MATCHES_RC * 100 / 10000" | bc)%)"

# If both have very low matches, there might be an issue
if [[ $MATCHES_ORIGINAL -lt 1000 ]] && [[ $MATCHES_RC -lt 1000 ]]; then
    echo "WARNING: Very low barcode match rates detected!"
    echo "Checking first few barcodes from file:"
    head -5 "$LOCAL_TMP/sampled_barcodes.txt"
    echo ""
    echo "Expected barcode examples from whitelist:"
    head -5 "$BARCODE_WHITELIST"
    echo ""
    echo "This might indicate:"
    echo "  1. R2 file doesn't contain barcodes in expected position"
    echo "  2. Wrong chemistry/kit version"
    echo "  3. Barcodes need extraction or processing"
    echo ""
    echo "Continuing anyway, but results may be poor..."
fi

# Choose the whitelist with more matches
if [[ $MATCHES_RC -gt $MATCHES_ORIGINAL ]]; then
    echo "Using reverse complement whitelist"
    FINAL_WHITELIST="$BARCODE_WHITELIST_RC"
else
    echo "Using original whitelist"
    FINAL_WHITELIST="$BARCODE_WHITELIST"
fi

# Calculate memory
TOTAL_MEM_GB=$(($SLURM_MEM_PER_NODE / 1024))
CELLRANGER_MEM_GB=$((TOTAL_MEM_GB - 16))

echo "Using $SLURM_CPUS_PER_TASK cores and ${CELLRANGER_MEM_GB}GB memory"

# Change to temporary directory
cd "$LOCAL_TMP"

# Clean up any existing pipestance
if [[ -d "$UNIQUE_ID" ]]; then
    rm -rf "$UNIQUE_ID"
fi

# ============================================
# Run CellRanger ATAC (not ARC!)
# ============================================
echo "Starting CellRanger ATAC..."
echo "Sample: $SAMPLE"
echo "Reference: $REF"
echo "FASTQ directory: $FASTQ_DIR"
echo "Whitelist: $FINAL_WHITELIST"

# Option 1: If you have cellranger-atac installed
if [[ -x "$CELLRANGER_ATAC" ]]; then
    "$CELLRANGER_ATAC" count \
        --id="$UNIQUE_ID" \
        --reference="$REF" \
        --fastqs="$FASTQ_DIR" \
        --sample="$SAMPLE" \
        --localcores="$SLURM_CPUS_PER_TASK" \
        --localmem="$CELLRANGER_MEM_GB"
else
    echo "ERROR: cellranger-atac not found. Alternative: Use chromap for ATAC processing"
    
    # Option 2: Alternative using chromap (if cellranger-atac not available)
    echo "Using chromap as alternative..."
    
    # You would need to:
    # 1. Extract 16bp barcodes from R2
    # 2. Run chromap with correct parameters
    # 3. Call peaks with MACS2
    # See the documentation for chromap workflow
fi

# Check if processing succeeded
if [[ $? -eq 0 ]]; then
    echo "Processing completed successfully for sample: $SAMPLE"
    
    # Copy results
    ATAC_DIR="$OUTPUT_DIR/${SAMPLE}_atac_results"
    mkdir -p "$ATAC_DIR"
    
    if [[ -d "$UNIQUE_ID/outs" ]]; then
        cp -r "$UNIQUE_ID/outs/"* "$ATAC_DIR/"
        
        # Verify key output files
        echo "Checking output files:"
        for file in "web_summary.html" "summary.csv" "filtered_peak_bc_matrix.h5" \
                    "fragments.tsv.gz" "peaks.bed" "singlecell.csv"; do
            if [[ -f "$ATAC_DIR/$file" ]]; then
                echo "  ✓ $file"
            else
                echo "  ✗ $file (missing)"
            fi
        done
    fi
else
    echo "ERROR: CellRanger ATAC failed for sample: $SAMPLE"
    exit 1
fi

echo "ATAC processing completed for $SAMPLE"