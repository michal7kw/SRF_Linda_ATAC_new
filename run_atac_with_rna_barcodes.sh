#!/bin/bash
#SBATCH --job-name=atac_rna_barcodes
#SBATCH --output=logs/atac_rna_barcodes_%A_%a.out
#SBATCH --error=logs/atac_rna_barcodes_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# Process ATAC data using barcodes from successful RNA processing
# This handles the case where RNA and ATAC were sequenced separately

set -euo pipefail

# Load Python module if needed
module load python/3.8 2>/dev/null || true

# Configuration
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/nestin"
RNA_RESULTS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_RNA/cellranger_final_count_data_trans"
OUTPUT_BASE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/rna_barcode_filtered"
SCRIPTS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data"

# Sample names
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing ATAC sample using RNA barcodes: $SAMPLE"

# Create output directories
OUTPUT_DIR="$OUTPUT_BASE/$SAMPLE"
mkdir -p "$OUTPUT_DIR" "logs"

# Step 1: Extract barcodes from RNA if not already done
BARCODE_DIR="$OUTPUT_BASE/barcode_whitelists"
mkdir -p "$BARCODE_DIR"

RNA_SAMPLE_NAME="${SAMPLE//-/_}"  # Convert format
BARCODE_FILE="$BARCODE_DIR/${SAMPLE}_rna_barcodes.txt"

if [[ ! -f "$BARCODE_FILE" ]]; then
    echo "Extracting barcodes from RNA results..."
    
    # Look for RNA matrix
    MATRIX_DIR=""
    for dir in \
        "$RNA_RESULTS_DIR/cellranger_counts_${RNA_SAMPLE_NAME}/outs/filtered_feature_bc_matrix" \
        "$RNA_RESULTS_DIR/${RNA_SAMPLE_NAME}_counts/filtered_feature_bc_matrix" \
        "$RNA_RESULTS_DIR/${RNA_SAMPLE_NAME}/outs/filtered_feature_bc_matrix"; do
        if [[ -d "$dir" ]]; then
            MATRIX_DIR="$dir"
            break
        fi
    done
    
    if [[ -z "$MATRIX_DIR" ]]; then
        echo "ERROR: Could not find RNA matrix for $SAMPLE"
        echo "Searched in $RNA_RESULTS_DIR"
        exit 1
    fi
    
    # Extract barcodes
    if [[ -f "$MATRIX_DIR/barcodes.tsv.gz" ]]; then
        zcat "$MATRIX_DIR/barcodes.tsv.gz" | cut -f1 -d'-' > "$BARCODE_FILE"
    elif [[ -f "$MATRIX_DIR/barcodes.tsv" ]]; then
        cut -f1 -d'-' "$MATRIX_DIR/barcodes.tsv" > "$BARCODE_FILE"
    else
        echo "ERROR: No barcodes file found in $MATRIX_DIR"
        exit 1
    fi
    
    NUM_BARCODES=$(wc -l < "$BARCODE_FILE")
    echo "Extracted $NUM_BARCODES valid cell barcodes from RNA"
else
    NUM_BARCODES=$(wc -l < "$BARCODE_FILE")
    echo "Using existing barcode file with $NUM_BARCODES barcodes"
fi

# Step 2: Analyze barcode positions in ATAC R2
echo "Analyzing barcode positions in ATAC data..."
python3 "$SCRIPTS_DIR/process_atac_custom_barcodes.py" \
    --sample "$SAMPLE" \
    --data-dir "$DATA_DIR" \
    --whitelist "$BARCODE_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --analyze-only

# Step 3: Filter ATAC reads using RNA barcodes
echo "Filtering ATAC reads with valid RNA barcodes..."
python3 "$SCRIPTS_DIR/process_atac_custom_barcodes.py" \
    --sample "$SAMPLE" \
    --data-dir "$DATA_DIR" \
    --whitelist "$BARCODE_FILE" \
    --output-dir "$OUTPUT_DIR"

# Step 4: Process filtered reads with standard ATAC pipeline
FILTERED_DIR="$OUTPUT_DIR"
FINAL_OUTPUT="$OUTPUT_BASE/${SAMPLE}_processed"
mkdir -p "$FINAL_OUTPUT"

# Check if we have filtered files
if [[ -f "$FILTERED_DIR/${SAMPLE}_filtered_R1_001.fastq.gz" ]]; then
    echo "Filtered FASTQ files created successfully"
    
    # Now we can try processing with CellRanger ATAC on the filtered data
    # Since we've pre-filtered to valid cells, it might work better
    
    # Create a simple barcode whitelist for CellRanger
    WHITELIST_10X="$FINAL_OUTPUT/barcodes_10x.txt"
    cat "$BARCODE_FILE" | awk '{print $1"-1"}' > "$WHITELIST_10X"
    
    echo "Summary:"
    echo "- Original ATAC data: $DATA_DIR"
    echo "- RNA barcodes used: $NUM_BARCODES cells"
    echo "- Filtered ATAC data: $FILTERED_DIR"
    echo "- Barcode whitelist: $WHITELIST_10X"
    
    # Optional: Run a simple peak calling pipeline
    echo ""
    echo "Next steps:"
    echo "1. Use filtered FASTQ files for alignment with BWA or bowtie2"
    echo "2. Call peaks with MACS2"
    echo "3. Create count matrix with ArchR or custom scripts"
    echo "4. Integrate with RNA data using Seurat or Signac"
    
    # Save processing info
    cat > "$FINAL_OUTPUT/processing_info.txt" <<EOF
Sample: $SAMPLE
Date: $(date)
RNA barcodes: $NUM_BARCODES
RNA source: $MATRIX_DIR
ATAC source: $DATA_DIR
Filtered output: $FILTERED_DIR
Processing script: $0
EOF
    
else
    echo "ERROR: Filtered files not created"
    exit 1
fi

echo "Processing complete for $SAMPLE"
echo "Filtered ATAC data ready for downstream analysis"