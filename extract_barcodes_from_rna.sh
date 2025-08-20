#!/bin/bash
#SBATCH --job-name=extract_rna_barcodes
#SBATCH --output=logs/extract_rna_barcodes_%A.out
#SBATCH --error=logs/extract_rna_barcodes_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=workq

# Extract valid cell barcodes from successful RNA processing
# These can be used as a whitelist for ATAC processing

set -euo pipefail

# Configuration - Point to successful RNA results
RNA_RESULTS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_RNA/cellranger_final_count_data"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Samples to process
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

echo "Extracting cell barcodes from successful RNA processing..."

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"
    
    # Find the RNA results directory (might have different naming)
    RNA_SAMPLE_NAME="${SAMPLE//-/_}"  # Convert R26-Nestin to R26_Nestin format
    
    # Look for the filtered matrix directory
    MATRIX_DIR=""
    if [[ -d "$RNA_RESULTS_DIR/cellranger_counts_${RNA_SAMPLE_NAME}/outs/filtered_feature_bc_matrix" ]]; then
        MATRIX_DIR="$RNA_RESULTS_DIR/cellranger_counts_${RNA_SAMPLE_NAME}/outs/filtered_feature_bc_matrix"
    elif [[ -d "$RNA_RESULTS_DIR/${RNA_SAMPLE_NAME}_counts/filtered_feature_bc_matrix" ]]; then
        MATRIX_DIR="$RNA_RESULTS_DIR/${RNA_SAMPLE_NAME}_counts/filtered_feature_bc_matrix"
    else
        echo "  WARNING: Could not find RNA matrix for $SAMPLE"
        echo "  Searched in:"
        echo "    - $RNA_RESULTS_DIR/cellranger_counts_${RNA_SAMPLE_NAME}/outs/filtered_feature_bc_matrix"
        echo "    - $RNA_RESULTS_DIR/${RNA_SAMPLE_NAME}_counts/filtered_feature_bc_matrix"
        continue
    fi
    
    # Extract barcodes from the matrix
    if [[ -f "$MATRIX_DIR/barcodes.tsv.gz" ]]; then
        echo "  Found barcodes file: $MATRIX_DIR/barcodes.tsv.gz"
        
        # Decompress and save barcodes
        OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE}_rna_barcodes.txt"
        zcat "$MATRIX_DIR/barcodes.tsv.gz" | cut -f1 -d'-' > "$OUTPUT_FILE"
        
        NUM_BARCODES=$(wc -l < "$OUTPUT_FILE")
        echo "  Extracted $NUM_BARCODES valid cell barcodes"
        echo "  Saved to: $OUTPUT_FILE"
        
        # Also create a 10x-formatted whitelist (with -1 suffix)
        OUTPUT_10X="$OUTPUT_DIR/${SAMPLE}_rna_barcodes_10x.txt"
        zcat "$MATRIX_DIR/barcodes.tsv.gz" > "$OUTPUT_10X"
        echo "  Also saved 10x format to: $OUTPUT_10X"
        
    elif [[ -f "$MATRIX_DIR/barcodes.tsv" ]]; then
        echo "  Found uncompressed barcodes file: $MATRIX_DIR/barcodes.tsv"
        
        # Copy barcodes
        OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE}_rna_barcodes.txt"
        cut -f1 -d'-' "$MATRIX_DIR/barcodes.tsv" > "$OUTPUT_FILE"
        
        NUM_BARCODES=$(wc -l < "$OUTPUT_FILE")
        echo "  Extracted $NUM_BARCODES valid cell barcodes"
        echo "  Saved to: $OUTPUT_FILE"
        
        # Also create a 10x-formatted whitelist
        OUTPUT_10X="$OUTPUT_DIR/${SAMPLE}_rna_barcodes_10x.txt"
        cp "$MATRIX_DIR/barcodes.tsv" "$OUTPUT_10X"
        echo "  Also saved 10x format to: $OUTPUT_10X"
    else
        echo "  ERROR: No barcodes file found in $MATRIX_DIR"
    fi
done

echo ""
echo "Barcode extraction complete!"
echo "Whitelists saved in: $OUTPUT_DIR"
echo ""
echo "These barcodes can be used for:"
echo "1. Filtering ATAC reads to only valid cells"
echo "2. Creating a custom barcode whitelist for CellRanger"
echo "3. Matching ATAC and RNA data from the same cell types"