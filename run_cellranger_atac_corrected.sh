#!/bin/bash
#SBATCH --job-name=cellranger_atac_corrected
#SBATCH --output=logs/cellranger_atac_corrected_%A_%a.out
#SBATCH --error=logs/cellranger_atac_corrected_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --partition=workq

# Final Corrected CellRanger ATAC Processing Script
# This script resolves the barcode issue by swapping the R1 and R2 FASTQ files,
# placing the barcode information where cellranger-atac expects it.

set -euo pipefail

# Configuration
CELLRANGER_ATAC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-atac-2.2.0/cellranger-atac"
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A"
SOURCE_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/nestin"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_atac_corrected_output"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp"

# Sample names
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get sample for this array job
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing ATAC sample with corrected read structure: $SAMPLE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "logs"

# Create a temporary directory for the renamed FASTQ files
RENAMED_FASTQ_DIR="$TMP_DIR/renamed_fastqs_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$RENAMED_FASTQ_DIR"

echo "Creating symbolic links with swapped R1 and R2 reads..."
# The key insight: R2 contains the barcodes, R1 contains genomic DNA.
# We swap them to match the standard 10x ATAC format.
ln -sf "$SOURCE_DATA_DIR/${SAMPLE}_R2_001.fastq.gz" "$RENAMED_FASTQ_DIR/${SAMPLE}_S1_L001_R1_001.fastq.gz" # R2 -> R1 (barcodes)
ln -sf "$SOURCE_DATA_DIR/${SAMPLE}_R1_001.fastq.gz" "$RENAMED_FASTQ_DIR/${SAMPLE}_S1_L001_R2_001.fastq.gz" # R1 -> R2 (genomic)
ln -sf "$SOURCE_DATA_DIR/${SAMPLE}_R3_001.fastq.gz" "$RENAMED_FASTQ_DIR/${SAMPLE}_S1_L001_R3_001.fastq.gz" # R3 is read 2
ln -sf "$SOURCE_DATA_DIR/${SAMPLE}_I1_001.fastq.gz" "$RENAMED_FASTQ_DIR/${SAMPLE}_S1_L001_I1_001.fastq.gz" # I1 is index

echo "Starting CellRanger ATAC..."
UNIQUE_ID="${SAMPLE}_atac_corrected_$(date +%Y%m%d%H%M%S)"
LOCAL_TMP="$TMP_DIR/cellranger_atac_corrected_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$LOCAL_TMP"
cd "$LOCAL_TMP"

"$CELLRANGER_ATAC" count \
    --id="$UNIQUE_ID" \
    --reference="$REF" \
    --fastqs="$RENAMED_FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores="$SLURM_CPUS_PER_TASK" \
    --localmem=112

if [[ $? -eq 0 ]]; then
    echo "CellRanger ATAC processing completed successfully for sample: $SAMPLE"
    RESULTS_DIR="$OUTPUT_DIR/${SAMPLE}_results"
    mkdir -p "$RESULTS_DIR"
    cp -r "$UNIQUE_ID/outs/"* "$RESULTS_DIR/"
    echo "Results successfully copied to $RESULTS_DIR"
else
    echo "ERROR: CellRanger ATAC failed for sample: $SAMPLE"
    exit 1
fi

echo "Processing completed for $SAMPLE"