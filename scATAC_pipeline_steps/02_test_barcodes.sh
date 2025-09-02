#!/bin/bash
#SBATCH --job-name=test_barcodes
#SBATCH --output=logs/test_barcodes_%A_%a.out
#SBATCH --error=logs/test_barcodes_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 2: Testing barcode quality for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create output directories
mkdir -p "$OUTPUT_DIR/qc"

EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "ERROR: Extracted barcode file not found: $EXTRACTED_BC_FILE"
    echo "Please run step 1 (01_extract_barcodes.sh) first"
    exit 1
fi

# Download and prepare whitelists
WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"
mkdir -p "$WHITELIST_DIR"

# Download whitelists if needed
if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1.txt" ]]; then
    echo "DEBUG: Downloading whitelist..."
    wget -q -O "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz" \
        https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to download whitelist"
        exit 1
    fi
    
    gunzip "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz"
    echo "DEBUG: Whitelist downloaded and extracted"
fi

echo "DEBUG: Checking whitelist files..."
echo "  Original whitelist: $(wc -l < "$WHITELIST_DIR/atac_737K-arc-v1.txt") barcodes"

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" ]]; then
    echo "DEBUG: Creating reverse complement whitelist..."
    cat "$WHITELIST_DIR/atac_737K-arc-v1.txt" | rev | tr 'ACGT' 'TGCA' > \
        "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
fi

echo "  RC whitelist: $(wc -l < "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt") barcodes"

# Test extracted barcodes
echo "DEBUG: Testing extracted barcodes against whitelists..."
echo "DEBUG: Extracting test barcodes..."

# Extract test barcodes (should now be exactly 16bp)
timeout 300 bash -c "zcat '$EXTRACTED_BC_FILE' | awk 'NR%4==2' | head -10000 > '$OUTPUT_DIR/qc/test_barcodes.txt'"
EXTRACT_EXIT_CODE=$?

if [[ $EXTRACT_EXIT_CODE -eq 124 ]]; then
    echo "ERROR: Barcode extraction timed out"
    exit 1
elif [[ $EXTRACT_EXIT_CODE -ne 0 ]]; then
    echo "ERROR: Barcode extraction failed"
    exit 1
fi

TEST_BC_COUNT=$(wc -l < "$OUTPUT_DIR/qc/test_barcodes.txt")
echo "DEBUG: Extracted $TEST_BC_COUNT test barcodes"

if [[ $TEST_BC_COUNT -eq 0 ]]; then
    echo "ERROR: No test barcodes extracted"
    exit 1
fi

# Check barcode length
SAMPLE_SEQ=$(head -1 "$OUTPUT_DIR/qc/test_barcodes.txt")
SEQ_LENGTH=${#SAMPLE_SEQ}
echo "DEBUG: Sample barcode: $SAMPLE_SEQ"
echo "DEBUG: Barcode length: $SEQ_LENGTH"

if [[ $SEQ_LENGTH -ne 16 ]]; then
    echo "ERROR: Expected 16bp barcodes, got ${SEQ_LENGTH}bp"
    exit 1
fi

# Analyze N content
N_COUNT=$(echo "$SAMPLE_SEQ" | grep -o 'N' | wc -l)
echo "DEBUG: Number of N's in sample barcode: $N_COUNT"

# Test both whitelist orientations with timeout and optimization
echo "DEBUG: Testing whitelist orientations..."
echo "DEBUG: Using optimized matching approach to avoid hangs..."

# Use a smaller test set for initial orientation testing
head -1000 "$OUTPUT_DIR/qc/test_barcodes.txt" > "$OUTPUT_DIR/qc/test_barcodes_small.txt"

echo "DEBUG: Testing original whitelist..."
MATCHES_ORIG=$(timeout 60 bash -c "grep -Fxf '$WHITELIST_DIR/atac_737K-arc-v1.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt' | wc -l" 2>/dev/null || echo "0")
if [[ $? -eq 124 ]]; then
    echo "WARNING: Original whitelist test timed out, using alternative method"
    MATCHES_ORIG=$(timeout 30 bash -c "awk 'NR==FNR{a[\$0]=1;next} \$0 in a{c++} END{print c+0}' '$WHITELIST_DIR/atac_737K-arc-v1.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt'" 2>/dev/null || echo "0")
fi

echo "DEBUG: Testing reverse complement whitelist..."
MATCHES_RC=$(timeout 60 bash -c "grep -Fxf '$WHITELIST_DIR/atac_737K-arc-v1_rc.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt' | wc -l" 2>/dev/null || echo "0")
if [[ $? -eq 124 ]]; then
    echo "WARNING: Reverse complement whitelist test timed out, using alternative method"
    MATCHES_RC=$(timeout 30 bash -c "awk 'NR==FNR{a[\$0]=1;next} \$0 in a{c++} END{print c+0}' '$WHITELIST_DIR/atac_737K-arc-v1_rc.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt'" 2>/dev/null || echo "0")
fi

SMALL_BC_COUNT=$(wc -l < "$OUTPUT_DIR/qc/test_barcodes_small.txt")
echo "Whitelist matching results (tested on $SMALL_BC_COUNT barcodes):"
echo "  Original: $MATCHES_ORIG / $SMALL_BC_COUNT ($(echo "scale=2; ${MATCHES_ORIG:-0} * 100 / $SMALL_BC_COUNT" | bc)%)"
echo "  Reverse complement: $MATCHES_RC / $SMALL_BC_COUNT ($(echo "scale=2; ${MATCHES_RC:-0} * 100 / $SMALL_BC_COUNT" | bc)%)"

# Choose the best matching whitelist
if [[ $MATCHES_RC -gt $MATCHES_ORIG ]]; then
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
    BEST_WHITELIST="reverse_complement"
    echo "Using reverse complement whitelist"
else
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1.txt"
    BEST_WHITELIST="original"
    echo "Using original whitelist"
fi

# Set barcode error threshold based on match rate
MAX_MATCHES=$((MATCHES_ORIG > MATCHES_RC ? MATCHES_ORIG : MATCHES_RC))
MATCH_RATE=$(echo "scale=2; ${MAX_MATCHES:-0} * 100 / $SMALL_BC_COUNT" | bc)

echo "DEBUG: Best match rate: ${MATCH_RATE}% (from $SMALL_BC_COUNT test barcodes)"

if [[ $(echo "$MATCH_RATE < 5.0" | bc) -eq 1 ]]; then
    echo "ERROR: Very low barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=3
    exit 1
elif [[ $(echo "$MATCH_RATE < 20.0" | bc) -eq 1 ]]; then
    echo "WARNING: Low barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=2
elif [[ $(echo "$MATCH_RATE < 50.0" | bc) -eq 1 ]]; then
    echo "INFO: Moderate barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=1
else
    echo "INFO: Good barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=1
fi

# Save results for next steps
cat > "$OUTPUT_DIR/qc/${SAMPLE}_barcode_test_results.txt" << EOF
SAMPLE=$SAMPLE
WHITELIST=$WHITELIST
BEST_WHITELIST=$BEST_WHITELIST
MATCH_RATE=$MATCH_RATE
BC_ERROR_THRESHOLD=$BC_ERROR_THRESHOLD
MATCHES_ORIG=$MATCHES_ORIG
MATCHES_RC=$MATCHES_RC
TEST_BC_COUNT=$SMALL_BC_COUNT
SAMPLE_BARCODE=$SAMPLE_SEQ
N_COUNT=$N_COUNT
EOF

echo "Results saved to: $OUTPUT_DIR/qc/${SAMPLE}_barcode_test_results.txt"

echo "========================================="
echo "Step 2 complete for $SAMPLE"
echo "Best whitelist: $BEST_WHITELIST (${MATCH_RATE}% match rate)"
echo "End time: $(date)"
echo "========================================="