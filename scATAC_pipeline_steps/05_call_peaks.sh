#!/bin/bash
#SBATCH --job-name=call_peaks
#SBATCH --output=logs/05_call_peaks_%a.out
#SBATCH --error=logs/05_call_peaks_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"
REF_GENOME="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"

echo "========================================="
echo "Step 5: Peak calling for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create output directories
mkdir -p "$OUTPUT_DIR/peaks" "$OUTPUT_DIR/qc"

# Check prerequisites
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    echo "Please run step 4 (04_chromap_alignment.sh) first"
    exit 1
fi

# Generate genome sizes
echo "DEBUG: Checking for genome sizes file..."
if [[ ! -f "$OUTPUT_DIR/qc/genome.chrom.sizes" ]]; then
    echo "DEBUG: Generating genome sizes file from reference index..."
    if [[ -f "${REF_GENOME}.fai" ]]; then
        cut -f1,2 "${REF_GENOME}.fai" > "$OUTPUT_DIR/qc/genome.chrom.sizes"
        echo "DEBUG: Created genome sizes file: $(wc -l < "$OUTPUT_DIR/qc/genome.chrom.sizes") chromosomes"
    else
        echo "DEBUG: Creating reference index first..."
        samtools faidx "$REF_GENOME"
        cut -f1,2 "${REF_GENOME}.fai" > "$OUTPUT_DIR/qc/genome.chrom.sizes"
        echo "DEBUG: Created genome sizes file: $(wc -l < "$OUTPUT_DIR/qc/genome.chrom.sizes") chromosomes"
    fi
else
    echo "DEBUG: Genome sizes file already exists: $(wc -l < "$OUTPUT_DIR/qc/genome.chrom.sizes") chromosomes"
fi

# Convert fragments to reads for peak calling
echo "DEBUG: Converting fragments to reads for peak calling..."

if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    # Check if bedClip is available
    # A more robust way to generate reads from fragments, avoiding complex awk scripts
    echo "DEBUG: Generating reads from fragment ends..."
    zcat "$FRAGMENTS_FILE" | \
        grep -v "^chrM" | \
        awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1; print $1, $3-1, $3}' | \
        sort -k1,1 -k2,2n | \
        gzip > "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"

    READ_COUNT=$(zcat "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" | wc -l)

    if [[ $? -eq 0 && $READ_COUNT -gt 0 ]]; then
        echo "DEBUG: Reads file created successfully"
        echo "DEBUG: Number of reads: $(printf "%'d" $READ_COUNT)"
    else
        echo "ERROR: Failed to create reads file or the file is empty"
        echo "DEBUG: Read count: $READ_COUNT"
        exit 1
    fi

    # Call peaks
    echo "DEBUG: Checking for MACS2..."
    if ! command -v macs2 &> /dev/null; then
        echo "ERROR: MACS2 not found. Please install MACS2 for peak calling"
        echo "DEBUG: You can install with: pip install MACS2"
        exit 1
    fi

    echo "DEBUG: MACS2 found: $(which macs2)"
    echo "DEBUG: MACS2 version: $(macs2 --version 2>&1 || echo 'Version check failed')"

    echo "DEBUG: Running MACS2 peak calling..."
    macs2 callpeak \
        -t "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" \
        -g mm -f BED -q 0.01 \
        --nomodel --shift -100 --extsize 200 \
        --keep-dup all \
        -B --SPMR \
        --outdir "$OUTPUT_DIR/peaks" \
        -n "${SAMPLE}"
        
    MACS2_EXIT_CODE=$?
    echo "DEBUG: MACS2 exit code: $MACS2_EXIT_CODE"
    
    if [[ $MACS2_EXIT_CODE -ne 0 ]]; then
        echo "ERROR: MACS2 failed with exit code $MACS2_EXIT_CODE"
        exit 1
    fi
else
    echo "DEBUG: Skipping MACS2 peak calling, peaks file already exists."
fi

if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    echo "ERROR: MACS2 completed but peaks file not found"
    exit 1
fi

# Peak statistics
PEAK_COUNT=$(wc -l < "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak")
echo "DEBUG: Number of peaks: $PEAK_COUNT"

# Create sorted peaks file for downstream analysis
echo "DEBUG: Creating sorted peaks file..."
cut -f 1-4 "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" | \
    sort -k1,1 -k2,2n > "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

echo "Peak calling results:"
echo "  Total peaks called: $(printf "%'d" $PEAK_COUNT)"
echo "  Peak file: $OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak"
echo "  Sorted peaks: $OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

# Additional MACS2 output files
echo "Additional MACS2 output files:"
for ext in summits.bed control_lambda.bdg treat_pileup.bdg; do
    MACS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_${ext}"
    if [[ -f "$MACS_FILE" ]]; then
        echo "  - $(basename "$MACS_FILE"): $(wc -l < "$MACS_FILE") lines"
    fi
done

# Save peak calling results
cat > "$OUTPUT_DIR/peaks/${SAMPLE}_peak_calling_results.txt" << EOF
SAMPLE=$SAMPLE
PEAK_COUNT=$PEAK_COUNT
FRAGMENTS_INPUT=$FRAGMENTS_FILE
PEAK_FILE=$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak
SORTED_PEAKS=$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed
PEAK_CALLING_COMPLETED_AT=$(date)
MACS2_PARAMETERS=-g mm -f BED -q 0.01 --nomodel --shift -100 --extsize 200 --keep-dup all -B --SPMR
EOF

echo "Results saved to: $OUTPUT_DIR/peaks/${SAMPLE}_peak_calling_results.txt"

echo "========================================="
echo "Step 5 complete for $SAMPLE"
echo "Called $(printf "%'d" $PEAK_COUNT) peaks"
echo "End time: $(date)"
echo "========================================="