#!/bin/bash
#SBATCH --job-name=cellranger_atac_nestin
#SBATCH --output=logs/cellranger_atac_%A_%a.out
#SBATCH --error=logs/cellranger_atac_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Set up cleanup trap
cleanup() {
    echo "Cleaning up processes..."
    killall -9 cellranger-atac 2>/dev/null
    killall -9 cellranger 2>/dev/null
    echo "Cleanup complete."
    exit
}

# Prevent hard link errors
export CELLRANGER_COPY_MODE=copy
export CELLRANGER_USE_HARDLINKS=false

# Trap signals
trap cleanup SIGINT SIGTERM EXIT

# Create logs directory if it doesn't exist
mkdir -p logs

# Path to Cell Ranger ATAC (for ATAC-only processing)
CELLRANGER_ATAC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-atac-2.2.0/cellranger-atac"

# Path to ATAC data directory
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/nestin"

# Path to ATAC reference genome (mouse GRCm39 - compatible with ATAC)
REF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A"

# Output directory
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/cellranger_atac_output"
mkdir -p $OUTPUT_DIR

# Sample names for Nestin samples
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Get the current sample based on array task ID
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing ATAC sample: $SAMPLE"

# Generate a timestamp for unique IDs
TIMESTAMP=$(date +%Y%m%d%H%M%S)

# Create a unique ID for this run
UNIQUE_ID="${SAMPLE}_atac_${TIMESTAMP}"

# Create a directory for fastq files
FASTQ_DIR="$OUTPUT_DIR/${SAMPLE}_fastq"
mkdir -p $FASTQ_DIR

# Create symbolic links for ATAC files
# ATAC-seq typically has I1 (index), R1 (cell barcode), R2 and R3 (genomic DNA)
echo "Creating symbolic links for ATAC files..."
ln -sf "$DATA_DIR/${SAMPLE}_I1_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_I1_001.fastq.gz"
ln -sf "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_R1_001.fastq.gz"
ln -sf "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_R2_001.fastq.gz"
ln -sf "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" "$FASTQ_DIR/${SAMPLE}_S1_L001_R3_001.fastq.gz"

# Verify files exist
echo "Verifying input files exist:"
for file in "$DATA_DIR/${SAMPLE}_I1_001.fastq.gz" \
           "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" \
           "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" \
           "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    else
        echo "  Found: $file"
    fi
done

# Get number of physical cores on the node for optimal performance
CORES=$SLURM_CPUS_PER_TASK
# Calculate memory in GB, leaving some overhead for the system
MEM_GB=$((${SLURM_MEM_PER_NODE} / 1024 - 16))

echo "Using $CORES cores and ${MEM_GB}GB memory for Cell Ranger ATAC"

# Add this function before running cellranger
ensure_directory() {
    local max_retries=5
    local count=0
    while [ $count -lt $max_retries ]; do
        mkdir -p "$1" && break
        echo "Failed to create directory $1, retrying in 10s..."
        sleep 10
        count=$((count+1))
    done
    if [ $count -eq $max_retries ]; then
        echo "Failed to create directory after $max_retries attempts"
        return 1
    fi
    return 0
}

# Then use it before each mkdir
ensure_directory $OUTPUT_DIR
ensure_directory $FASTQ_DIR

# Create a local temporary directory
LOCAL_TMP="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp/cellranger_atac_${SAMPLE}_${SLURM_ARRAY_TASK_ID}"
mkdir -p $LOCAL_TMP

# Clean up any existing pipestance directory that might conflict
if [ -d "${UNIQUE_ID}" ]; then
    echo "Removing existing pipestance directory: ${UNIQUE_ID}"
    rm -rf "${UNIQUE_ID}"
fi

# Also check in the temp directory
if [ -d "${LOCAL_TMP}/${UNIQUE_ID}" ]; then
    echo "Removing existing pipestance directory in temp location: ${LOCAL_TMP}/${UNIQUE_ID}"
    rm -rf "${LOCAL_TMP}/${UNIQUE_ID}"
fi

# Check if Cell Ranger ATAC executable exists
if [ ! -f "$CELLRANGER_ATAC" ]; then
    echo "ERROR: Cell Ranger ATAC not found at: $CELLRANGER_ATAC"
    echo "Please verify the path or install Cell Ranger ATAC"
    exit 1
fi

# Check if reference genome exists
if [ ! -d "$REF" ]; then
    echo "ERROR: Reference genome not found at: $REF"
    echo "Please verify the path or download the ATAC reference genome"
    exit 1
fi

echo "Starting Cell Ranger ATAC processing..."
echo "Sample: $SAMPLE"
echo "Reference: $REF"
echo "FASTQ directory: $FASTQ_DIR"
echo "Output directory: $LOCAL_TMP"

# Run Cell Ranger ATAC count command - no libraries CSV needed for ATAC-only
cd $LOCAL_TMP
$CELLRANGER_ATAC count \
    --id=${UNIQUE_ID} \
    --reference=$REF \
    --fastqs=$FASTQ_DIR \
    --sample=$SAMPLE \
    --localcores=$CORES \
    --localmem=$MEM_GB

# Check if Cell Ranger ATAC completed successfully
if [ $? -ne 0 ]; then
    echo "ERROR: Cell Ranger ATAC failed for sample: $SAMPLE"
    echo "Check logs for more details"
    exit 1
fi

# Only proceed with copying if the output directory exists
if [ -d "$LOCAL_TMP/${UNIQUE_ID}/outs" ]; then
    echo "Creating a separate directory for ATAC results"
    ATAC_DIR="$OUTPUT_DIR/${SAMPLE}_atac_results"
    ensure_directory $ATAC_DIR

    # Copy essential ATAC-seq output files
    echo "Copying ATAC-seq output files..."
    
    # Peak-barcode matrices (filtered and raw)
    if [ -d "$LOCAL_TMP/${UNIQUE_ID}/outs/filtered_peak_bc_matrix" ]; then
        rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/filtered_peak_bc_matrix/ $ATAC_DIR/filtered_peak_bc_matrix/
    fi
    if [ -d "$LOCAL_TMP/${UNIQUE_ID}/outs/raw_peak_bc_matrix" ]; then
        rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/raw_peak_bc_matrix/ $ATAC_DIR/raw_peak_bc_matrix/
    fi
    
    # Summary and metrics files
    [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/singlecell.csv" ] && rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/singlecell.csv $ATAC_DIR/
    [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/summary.csv" ] && rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/summary.csv $ATAC_DIR/
    [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/web_summary.html" ] && rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/web_summary.html $ATAC_DIR/
    
    # Peak annotation files
    [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/peaks.bed" ] && rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/peaks.bed $ATAC_DIR/
    [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/peak_annotation.tsv" ] && rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/peak_annotation.tsv $ATAC_DIR/
    
    # Fragments file (important for downstream analysis)
    if [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/fragments.tsv.gz" ]; then
        rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/fragments.tsv.gz* $ATAC_DIR/
    fi
    
    # Position-sorted BAM file
    if [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/possorted_bam.bam" ]; then
        rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/possorted_bam.bam* $ATAC_DIR/
    fi
    
    # Cloupe file for Loupe Browser visualization
    if [ -f "$LOCAL_TMP/${UNIQUE_ID}/outs/cloupe.cloupe" ]; then
        rsync -av $LOCAL_TMP/${UNIQUE_ID}/outs/cloupe.cloupe $ATAC_DIR/
    fi

    # Copy the full output directory for complete results
    rsync -av $LOCAL_TMP/${UNIQUE_ID}/ $OUTPUT_DIR/${UNIQUE_ID}/
    
    echo "Successfully copied ATAC-seq results to: $ATAC_DIR"
    echo "Full output available at: $OUTPUT_DIR/${UNIQUE_ID}"
else
    echo "ERROR: Cell Ranger ATAC output directory not found at $LOCAL_TMP/${UNIQUE_ID}/outs"
    echo "Cell Ranger ATAC likely failed. Check logs for details."
    exit 1
fi

# Clean up temporary directory to save space
echo "Cleaning up temporary files..."
rm -rf $LOCAL_TMP

echo "Cell Ranger ATAC processing complete for sample: $SAMPLE"
echo "ATAC results available at: $ATAC_DIR"
echo "Full output available at: $OUTPUT_DIR/${UNIQUE_ID}"