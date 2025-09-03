#!/bin/bash
#SBATCH --job-name=rerun_macs2
#SBATCH --output=logs/05a_rerun_macs2_%a.out
#SBATCH --error=logs/05a_rerun_macs2_%a.err
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
echo "Step 5a: Rerunning Peak calling for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Check prerequisites
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
if [[ ! -f "$READS_FILE" ]]; then
    echo "ERROR: Reads file not found: $READS_FILE"
    echo "Please run step 5 (05_call_peaks.sh) first to generate the reads file."
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
    -t "$READS_FILE" \
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
READS_INPUT=$READS_FILE
PEAK_FILE=$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak
SORTED_PEAKS=$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed
PEAK_CALLING_COMPLETED_AT=$(date)
MACS2_PARAMETERS=-g mm -f BED -q 0.01 --nomodel --shift -100 --extsize 200 --keep-dup all -B --SPMR
EOF

echo "Results saved to: $OUTPUT_DIR/peaks/${SAMPLE}_peak_calling_results.txt"

echo "========================================="
echo "Step 5a complete for $SAMPLE"
echo "Called $(printf "%'d" $PEAK_COUNT) peaks"
echo "End time: $(date)"
echo "========================================="