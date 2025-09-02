#!/bin/bash
#SBATCH --job-name=chromap_align
#SBATCH --output=logs/04_chromap_align_%a.out
#SBATCH --error=logs/04_chromap_align_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
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
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"
REF_GENOME="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"
CHROMAP_INDEX="$OUTPUT_DIR/chromap_index/mm10.index"
THREADS=16

echo "========================================="
echo "Step 4: Chromap alignment for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create output directories
mkdir -p "$OUTPUT_DIR/fragments" "$OUTPUT_DIR/logs"

EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

# Check prerequisites
if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "ERROR: Extracted barcode file not found: $EXTRACTED_BC_FILE"
    echo "Please run steps 1-3 first"
    exit 1
fi

# Load barcode test results
BARCODE_RESULTS="$OUTPUT_DIR/qc/${SAMPLE}_barcode_test_results.txt"
if [[ ! -f "$BARCODE_RESULTS" ]]; then
    echo "ERROR: Barcode test results not found: $BARCODE_RESULTS"
    echo "Please run step 2 (02_test_barcodes.sh) first"
    exit 1
fi

# Source the results
source "$BARCODE_RESULTS"

echo "DEBUG: Using barcode test results:"
echo "  Best whitelist: $BEST_WHITELIST"
echo "  Match rate: $MATCH_RATE%"
echo "  Error threshold: $BC_ERROR_THRESHOLD"

# Check chromap availability
echo "DEBUG: Checking chromap availability..."
if ! command -v chromap &> /dev/null; then
    echo "DEBUG: chromap not found in PATH, trying to load module..."
    
    # Try multiple module loading approaches
    if module load chromap 2>/dev/null; then
        echo "DEBUG: Successfully loaded chromap module"
    elif module load Chromap 2>/dev/null; then
        echo "DEBUG: Successfully loaded Chromap module (capitalized)"
    elif module load bioinformatics/chromap 2>/dev/null; then
        echo "DEBUG: Successfully loaded bioinformatics/chromap module"
    else
        echo "ERROR: chromap not found and module load failed"
        echo "DEBUG: Available modules containing 'chromap':"
        module avail chromap 2>&1 || echo "No chromap modules found"
        
        # Try conda installation
        if command -v conda &> /dev/null; then
            echo "DEBUG: Installing chromap via conda..."
            conda install -y -c bioconda chromap || {
                echo "ERROR: Failed to install chromap via conda"
                echo "Please install chromap manually from: https://github.com/haowenz/chromap"
                exit 1
            }
        else
            echo "ERROR: Neither chromap nor conda available"
            echo "Please install chromap from: https://github.com/haowenz/chromap"
            exit 1
        fi
    fi
    
    # Check again after module load/installation
    if ! command -v chromap &> /dev/null; then
        echo "ERROR: chromap still not available after module load/installation"
        exit 1
    fi
fi

echo "DEBUG: chromap found: $(which chromap)"
echo "DEBUG: chromap version: $(chromap --version 2>&1 | head -1 || echo 'Version check failed')"

# Check if chromap index exists, if not create it
if [[ ! -f "$CHROMAP_INDEX" ]]; then
    echo "DEBUG: Chromap index not found, creating index..."
    mkdir -p "$(dirname "$CHROMAP_INDEX")"
    
    echo "DEBUG: Building chromap index from reference: $REF_GENOME"
    chromap -i -r "$REF_GENOME" -o "$CHROMAP_INDEX"
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to build chromap index"
        exit 1
    fi
    echo "DEBUG: Chromap index created successfully: $CHROMAP_INDEX"
else
    echo "DEBUG: Using existing chromap index: $CHROMAP_INDEX"
fi

# Run chromap alignment
if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" ]]; then
    echo "DEBUG: Running chromap with the following parameters:"
    echo "  Reference: $REF_GENOME"
    echo "  Index: $CHROMAP_INDEX"
    echo "  R1: $DATA_DIR/${SAMPLE}_R1_001.fastq.gz"
    echo "  R2: $DATA_DIR/${SAMPLE}_R3_001.fastq.gz"
    echo "  Barcode: $EXTRACTED_BC_FILE (rightmost 16bp from R2)"
    echo "  Whitelist: $WHITELIST"
    echo "  Barcode error threshold: $BC_ERROR_THRESHOLD"
    echo "  Expected match rate: ${MATCH_RATE}%"

    chromap --preset atac \
        -x "$CHROMAP_INDEX" \
        -r "$REF_GENOME" \
        -1 "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" \
        -2 "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" \
        -b "$EXTRACTED_BC_FILE" \
        --barcode-whitelist "$WHITELIST" \
        --bc-error-threshold $BC_ERROR_THRESHOLD \
        -e 5 \
        -t $THREADS \
        --low-mem \
        -o "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv" 2>&1 | tee "$OUTPUT_DIR/logs/${SAMPLE}_chromap.log"
else
    echo "DEBUG: Skipping chromap alignment, fragments file already exists."
fi

# Check if chromap succeeded
if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv" && ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" ]]; then
    echo "ERROR: Chromap failed to generate fragments file"
    echo "DEBUG: Check chromap log: $OUTPUT_DIR/logs/${SAMPLE}_chromap.log"
    exit 1
fi

# Compress and index fragments
echo "DEBUG: Processing fragments..."
if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz.tbi" ]]; then
    bgzip -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv"
    tabix -f -s 1 -b 2 -e 3 -p bed "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
    echo "DEBUG: Fragments compressed and indexed"
else
    echo "DEBUG: Skipping fragment indexing, index file already exists."
fi

# Generate quick statistics
TOTAL_FRAGS=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | wc -l)
UNIQUE_BC=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | cut -f4 | sort -u | wc -l)
FRAGS_IN_CELLS=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | cut -f4 | grep -Ff "$WHITELIST" | wc -l)

echo "Fragment statistics:"
echo "  Total fragments: $TOTAL_FRAGS"
echo "  Unique barcodes: $UNIQUE_BC"
echo "  Fragments in cells: $FRAGS_IN_CELLS"
echo "  Fraction in cells: $(echo "scale=4; $FRAGS_IN_CELLS / $TOTAL_FRAGS" | bc)"

# Save alignment results
cat > "$OUTPUT_DIR/logs/${SAMPLE}_alignment_results.txt" << EOF
SAMPLE=$SAMPLE
TOTAL_FRAGS=$TOTAL_FRAGS
UNIQUE_BC=$UNIQUE_BC
FRAGS_IN_CELLS=$FRAGS_IN_CELLS
FRACTION_IN_CELLS=$(echo "scale=4; $FRAGS_IN_CELLS / $TOTAL_FRAGS" | bc)
WHITELIST_USED=$WHITELIST
MATCH_RATE_USED=$MATCH_RATE
BC_ERROR_THRESHOLD_USED=$BC_ERROR_THRESHOLD
ALIGNMENT_COMPLETED_AT=$(date)
EOF

echo "Results saved to: $OUTPUT_DIR/logs/${SAMPLE}_alignment_results.txt"

echo "========================================="
echo "Step 4 complete for $SAMPLE"
echo "Generated $(printf "%'d" $TOTAL_FRAGS) fragments from $(printf "%'d" $UNIQUE_BC) unique barcodes"
echo "End time: $(date)"
echo "========================================="