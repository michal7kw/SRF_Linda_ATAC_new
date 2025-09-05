#!/bin/bash
#SBATCH --job-name=peak_cell_matrix
#SBATCH --output=logs/06_peak_cell_matrix_%a.out
#SBATCH --error=logs/06_peak_cell_matrix_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

set -euo pipefail

# Helper function to zcat or cat a file based on gzip status
zcat_if_gzipped() {
    local file="$1"
    if file "$file" | grep -q "gzip compressed"; then
        zcat "$file"
    else
        cat "$file"
    fi
}
export -f zcat_if_gzipped # Export function for use in subshells/pipes

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 6: Peak-by-cell matrix for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Check prerequisites
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" # This file is expected to be gzipped
GENOME_SIZES="$OUTPUT_DIR/qc/genome.chrom.sizes"

if [[ ! -f "$PEAKS_FILE" ]]; then
    echo "ERROR: Sorted peaks file not found: $PEAKS_FILE"
    echo "Please run step 5 (05_call_peaks.sh) first"
    exit 1
fi

if [[ ! -f "$READS_FILE" ]]; then
    echo "ERROR: Reads file not found: $READS_FILE"
    echo "Please run step 5 (05_call_peaks.sh) first"
    exit 1
fi

if [[ ! -f "$GENOME_SIZES" ]]; then
    echo "ERROR: Genome sizes file not found: $GENOME_SIZES"
    exit 1
fi

# Generate peak-by-cell matrix
PEAK_OV_FILE="$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz"

# Force regeneration of PEAK_OV_FILE for debugging
rm -f "$PEAK_OV_FILE"

# Find reads in peaks per cell
echo "DEBUG: Checking for bedtools..."
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found. Please install bedtools"
    echo "DEBUG: You can install with: conda install -c bioconda bedtools"
    exit 1
fi

echo "DEBUG: bedtools found: $(which bedtools)"
echo "DEBUG: bedtools version: $(bedtools --version)"

echo "DEBUG: Finding reads overlapping with peaks per cell..."
echo "DEBUG: This may take several minutes for large datasets..."

echo "DEBUG: Inspecting PEAKS_FILE ($PEAKS_FILE) head:"
zcat_if_gzipped "$PEAKS_FILE" | head -n 5
echo "DEBUG: Inspecting READS_FILE ($READS_FILE) head:"
(zcat_if_gzipped "$READS_FILE" || true) | head -n 5

# The -sorted flag requires that both files are sorted by chrom, then start.
# The reads are piped from a sort command, and the peaks file was sorted in step 5.
bedtools intersect \
    -a "$PEAKS_FILE" \
    -b <(zcat_if_gzipped "$READS_FILE" | sort -k1,1 -k2,2n) \
    -wo -sorted -g "$GENOME_SIZES" | \
    tee "$OUTPUT_DIR/debug_${SAMPLE}_bedtools_intersect_raw.txt" | \
    head -n 10 > "$OUTPUT_DIR/debug_${SAMPLE}_bedtools_intersect_head.txt"

# Check if bedtools intersect produced any output
if [[ ! -s "$OUTPUT_DIR/debug_${SAMPLE}_bedtools_intersect_raw.txt" ]]; then
    echo "WARNING: bedtools intersect produced no overlaps for $SAMPLE."
    echo "========================================="
    echo "Step 6 complete for $SAMPLE (No data)"
    echo "End time: $(date)"
    echo "========================================="
    exit 0
fi

# Continue processing if overlaps were found
cat "$OUTPUT_DIR/debug_${SAMPLE}_bedtools_intersect_raw.txt" | \
    awk -v OFS='\t' '{print $4, $8}' | \
    tee "$OUTPUT_DIR/debug_${SAMPLE}_awk_output_raw.txt" | \
    head -n 10 > "$OUTPUT_DIR/debug_${SAMPLE}_awk_output_head.txt" | \
    sort -k1,1 -k2,2 | \
    uniq -c | \
    awk -v OFS='\t' '{print $2, $3, $1}' | \
    gzip > "$PEAK_OV_FILE"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: Peak-read overlap completed successfully"
    echo "DEBUG: Peak-read overlap file size: $(stat -c%s "$PEAK_OV_FILE") bytes"
    echo "DEBUG: Inspecting PEAK_OV_FILE ($PEAK_OV_FILE) head:"
    zcat_if_gzipped "$PEAK_OV_FILE" | head -n 10
else
    echo "ERROR: Failed to generate peak-read overlap"
    exit 1
fi

# Create matrix files in 10X format
echo "DEBUG: Creating 10X-compatible matrix files..."
MATRIX_DIR="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
mkdir -p "$MATRIX_DIR"

# 1. features.tsv (peaks)
# Use the 4th column (name) from the sorted peaks file as the feature ID.
# The second column (gene name) is the same for scATAC-seq.
# The third column is the type.
echo "DEBUG: Creating features.tsv file..."
cut -f 4 "$PEAKS_FILE" | awk -v OFS='\t' '{print $1, $1, "Peaks"}' > "$MATRIX_DIR/features.tsv"

# 2. barcodes.tsv
# Extract all unique barcodes from the overlap file.
echo "DEBUG: Creating barcodes.tsv file..."
zcat_if_gzipped "$PEAK_OV_FILE" | cut -f 2 | sort -u > "$MATRIX_DIR/barcodes.tsv"

BARCODE_COUNT=$(wc -l < "$MATRIX_DIR/barcodes.tsv")
if [[ $BARCODE_COUNT -eq 0 ]]; then
    echo "WARNING: No barcodes with overlapping reads found. Skipping matrix generation."
    echo "========================================="
    echo "Step 6 complete for $SAMPLE (No data)"
    echo "End time: $(date)"
    echo "========================================="
    exit 0
fi

# 3. matrix.mtx
# Convert the triplet file (peak_id, barcode_id, count) to an indexed MTX file.

echo "DEBUG: Creating indexed matrix.mtx file..."

# Create associative arrays to map IDs to 1-based indices.
# Then process the triplets file to print out the indexed matrix.
awk -v peak_file="$MATRIX_DIR/features.tsv" -v barcode_file="$MATRIX_DIR/barcodes.tsv" '
    BEGIN {OFS="\t"}
    # Read features into map
    FILENAME==peak_file {peak_map[$1] = FNR; next}
    # Read barcodes into map
    FILENAME==barcode_file {barcode_map[$1] = FNR; next}
    # Process triplets and print indexed matrix
    {
        peak_id = $1
        barcode_id = $2
        count = $3
        peak_idx = peak_map[peak_id]
        barcode_idx = barcode_map[barcode_id]
        if (peak_idx > 0 && barcode_idx > 0) {
            print peak_idx, barcode_idx, count
        }
    }' "$MATRIX_DIR/features.tsv" "$MATRIX_DIR/barcodes.tsv" <(zcat_if_gzipped "$PEAK_OV_FILE") > "$MATRIX_DIR/matrix.mtx.tmp"


TOTAL_PEAKS=$(wc -l < "$MATRIX_DIR/features.tsv")
TOTAL_BARCODES=$(wc -l < "$MATRIX_DIR/barcodes.tsv")
TOTAL_ENTRIES=$(wc -l < "$MATRIX_DIR/matrix.mtx.tmp")

# Create proper MTX format
echo "DEBUG: Creating MTX format matrix..."
{
    echo "%%MatrixMarket matrix coordinate integer general"
    echo "% Generated by scATAC pipeline step 6"
    echo "$TOTAL_PEAKS $TOTAL_BARCODES $TOTAL_ENTRIES"
    cat "$MATRIX_DIR/matrix.mtx.tmp"
} > "$MATRIX_DIR/matrix.mtx"

# Cleanup temporary file
rm -f "$MATRIX_DIR/matrix.mtx.tmp"

# Compress matrix files
echo "DEBUG: Compressing matrix files..."
gzip -f "$MATRIX_DIR/matrix.mtx"
gzip -f "$MATRIX_DIR/features.tsv"
gzip -f "$MATRIX_DIR/barcodes.tsv"

# Final stats and reporting
SPARSITY="N/A"
SPARSITY_VALUE="N/A"
if [[ $TOTAL_PEAKS -gt 0 && $TOTAL_BARCODES -gt 0 && $TOTAL_ENTRIES -gt 0 ]]; then
    SPARSITY_VALUE=$(echo "scale=4; (1 - $TOTAL_ENTRIES / ($TOTAL_PEAKS * $TOTAL_BARCODES)) * 100" | bc)
    SPARSITY="${SPARSITY_VALUE}%%"
fi

echo "Matrix statistics:"
echo "  Total peaks: $(printf "%d" $TOTAL_PEAKS)"
echo "  Total barcodes: $(printf "%d" $TOTAL_BARCODES)"
echo "  Matrix entries: $(printf "%d" $TOTAL_ENTRIES)"
echo "  Sparsity: ${SPARSITY}"

echo "Output files created:"
echo "  - $MATRIX_DIR/features.tsv.gz"
echo "  - $MATRIX_DIR/barcodes.tsv.gz"
echo "  - $MATRIX_DIR/matrix.mtx.gz"
echo "  - $PEAK_OV_FILE"

# Save matrix results
cat > "$MATRIX_DIR/matrix_info.txt" << EOF
SAMPLE=$SAMPLE
TOTAL_PEAKS=$TOTAL_PEAKS
TOTAL_BARCODES=$TOTAL_BARCODES
TOTAL_ENTRIES=$TOTAL_ENTRIES
SPARSITY="$SPARSITY"
SPARSITY_VALUE="$SPARSITY_VALUE"
MATRIX_CREATED_AT=$(date)
INPUT_PEAKS_FILE=$PEAKS_FILE
INPUT_READS_FILE=$READS_FILE
EOF

echo "Matrix info saved to: $MATRIX_DIR/matrix_info.txt"

echo "========================================="
echo "Step 6 complete for $SAMPLE"
echo "Created peak-by-cell matrix: $(printf "%d" $TOTAL_PEAKS) peaks × $(printf "%d" $TOTAL_BARCODES) cells"
echo "End time: $(date)"
echo "========================================="
