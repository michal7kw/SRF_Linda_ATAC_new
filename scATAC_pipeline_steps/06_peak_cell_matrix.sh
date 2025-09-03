#!/bin/bash
#SBATCH --job-name=peak_cell_matrix
#SBATCH --output=logs/06_peak_cell_matrix_%a.out
#SBATCH --error=logs/06_peak_cell_matrix_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

set -euo pipefail

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
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"

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

# Generate peak-by-cell matrix
if [[ ! -f "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" ]]; then
    
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
    
    bedtools intersect \
        -a "$PEAKS_FILE" \
        -b <(zcat "$READS_FILE" | sort -k1,1V -k2,2n) \
        -wo -sorted -g "$OUTPUT_DIR/qc/genome.chrom.sizes" | \
        sort -k8,8 | \
        bedtools groupby -g 8 -c 4 -o freqdesc | \
        gzip > "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz"

    if [[ $? -eq 0 ]]; then
        echo "DEBUG: Peak-read overlap completed successfully"
        echo "DEBUG: Peak-read overlap file size: $(stat -c%s "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz") bytes"
        echo "DEBUG: Number of cells with peaks: $(zcat "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" | wc -l)"
    else
        echo "ERROR: Failed to generate peak-read overlap"
        exit 1
    fi
else
    echo "DEBUG: Skipping peak-by-cell matrix generation, file already exists."
fi

# Create matrix files in 10X format
echo "DEBUG: Creating 10X-compatible matrix files..."
mkdir -p "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"

# Peaks bed file (features.tsv equivalent)
echo "DEBUG: Creating peaks.bed file..."
cut -f 1-3 "$PEAKS_FILE" > "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: peaks.bed created: $(wc -l < "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed") peaks"
else
    echo "ERROR: Failed to create peaks.bed"
    exit 1
fi

# Barcodes file
echo "DEBUG: Creating barcodes.tsv file..."
BARCODE_COUNT=$(zcat "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" | cut -f 1 | tee "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/barcodes.tsv" | wc -l)

if [[ $? -eq 0 ]]; then
    echo "DEBUG: barcodes.tsv created: $BARCODE_COUNT barcodes"
else
    echo "ERROR: Failed to create barcodes.tsv"
    exit 1
fi

if [[ $BARCODE_COUNT -eq 0 ]]; then
    echo "WARNING: No barcodes found. Skipping matrix generation."
    echo "========================================="
    echo "Step 6 complete for $SAMPLE (No data)"
    echo "End time: $(date)"
    echo "========================================="
    exit 0
fi

# Create features.tsv file (alternative to peaks.bed for 10X compatibility)
echo "DEBUG: Creating features.tsv file for 10X compatibility..."
awk 'BEGIN{OFS="\t"} {print $1":"$2"-"$3, $1":"$2"-"$3, "Peaks"}' \
    "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed" > \
    "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/features.tsv"

# Generate basic count matrix (simplified version)
echo "DEBUG: Creating basic peak-barcode count matrix..."
zcat "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" | \
    awk 'BEGIN{OFS="\t"} {
        split($2, peaks, ","); 
        split($3, counts, ","); 
        for(i=1; i<=length(peaks); i++) {
            if(counts[i] > 0) print $1, peaks[i], counts[i]
        }
    }' > "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix.mtx.tmp"

# Count unique peaks and barcodes for matrix header
TOTAL_PEAKS=$(wc -l < "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed")
TOTAL_BARCODES=$(wc -l < "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/barcodes.tsv")
TOTAL_ENTRIES=$(wc -l < "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix.mtx.tmp")

# Create proper MTX format
echo "DEBUG: Creating MTX format matrix..."
{
    echo "%%MatrixMarket matrix coordinate integer general"
    echo "% Generated by scATAC pipeline step 6"
    echo "$TOTAL_PEAKS $TOTAL_BARCODES $TOTAL_ENTRIES"
    cat "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix.mtx.tmp"
} > "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix.mtx"

# Cleanup temporary file
rm -f "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix.mtx.tmp"

# Compress matrix files
echo "DEBUG: Compressing matrix files..."
gzip -f "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix.mtx"
gzip -f "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/features.tsv"  
gzip -f "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/barcodes.tsv"

echo "Matrix statistics:"
echo "  Total peaks: $(printf "%'d" $TOTAL_PEAKS)"
echo "  Total barcodes: $(printf "%'d" $TOTAL_BARCODES)"
echo "  Matrix entries: $(printf "%'d" $TOTAL_ENTRIES)"
if [[ $TOTAL_PEAKS -gt 0 && $TOTAL_BARCODES -gt 0 ]]; then
    SPARSITY=$(echo "scale=4; (1 - $TOTAL_ENTRIES / ($TOTAL_PEAKS * $TOTAL_BARCODES)) * 100" | bc)
    echo "  Sparsity: ${SPARSITY}%"
else
    echo "  Sparsity: N/A (zero peaks or barcodes)"
fi

echo "Output files created:"
echo "  - $OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed"
echo "  - $OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/features.tsv.gz"
echo "  - $OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/barcodes.tsv.gz" 
echo "  - $OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix.mtx.gz"
echo "  - $OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz"

# Save matrix results
cat > "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix_info.txt" << EOF
SAMPLE=$SAMPLE
TOTAL_PEAKS=$TOTAL_PEAKS
TOTAL_BARCODES=$TOTAL_BARCODES
TOTAL_ENTRIES=$TOTAL_ENTRIES
SPARSITY="N/A"
if [[ $TOTAL_PEAKS -gt 0 && $TOTAL_BARCODES -gt 0 ]]; then
    SPARSITY_VALUE=$(echo "scale=4; (1 - $TOTAL_ENTRIES / ($TOTAL_PEAKS * $TOTAL_BARCODES)) * 100" | bc)
    SPARSITY="${SPARSITY_VALUE}%"
fi
MATRIX_CREATED_AT=$(date)
INPUT_PEAKS_FILE=$PEAKS_FILE
INPUT_READS_FILE=$READS_FILE
EOF

echo "Matrix info saved to: $OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/matrix_info.txt"

echo "========================================="
echo "Step 6 complete for $SAMPLE"
echo "Created peak-by-cell matrix: $(printf "%'d" $TOTAL_PEAKS) peaks × $(printf "%'d" $TOTAL_BARCODES) cells"
echo "End time: $(date)"
echo "========================================="