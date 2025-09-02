#!/bin/bash
#SBATCH --job-name=create_bigwig
#SBATCH --output=logs/07_create_bigwig_%a.out
#SBATCH --error=logs/07_create_bigwig_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 7: Creating BigWig files for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create bigwig output directory
mkdir -p "$OUTPUT_DIR/bigwig"

# Check prerequisites
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
CHROM_SIZES="$OUTPUT_DIR/mm10.chrom.sizes"

if [[ ! -f "$READS_FILE" ]]; then
    echo "ERROR: Reads file not found: $READS_FILE"
    echo "Please run chromap alignment first"
    exit 1
fi

if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "ERROR: Chromosome sizes file not found: $CHROM_SIZES"
    echo "Please ensure mm10.chrom.sizes exists in output directory"
    exit 1
fi

# Check for required tools
echo "DEBUG: Checking for required tools..."
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found. Please install bedtools"
    echo "DEBUG: You can install with: conda install -c bioconda bedtools"
    exit 1
fi

if ! command -v bedGraphToBigWig &> /dev/null; then
    echo "ERROR: bedGraphToBigWig not found. Please install UCSC tools"
    echo "DEBUG: You can install with: conda install -c bioconda ucsc-bedgraphtobigwig"
    exit 1
fi

echo "DEBUG: bedtools found: $(which bedtools)"
echo "DEBUG: bedGraphToBigWig found: $(which bedGraphToBigWig)"

# Create coverage bedGraph from reads
BEDGRAPH_FILE="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage.bedGraph"
BIGWIG_FILE="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage.bw"

if [[ ! -f "$BIGWIG_FILE" ]]; then
    echo "DEBUG: Creating coverage bedGraph..."
    echo "DEBUG: This may take several minutes for large datasets..."
    
    # Create bedGraph with coverage information
    # Using bedtools genomecov to create coverage track from BED file
    zcat "$READS_FILE" | \
        bedtools genomecov -i stdin -g "$CHROM_SIZES" -bg > "$BEDGRAPH_FILE"
    
    if [[ $? -eq 0 ]]; then
        echo "DEBUG: BedGraph created successfully"
        echo "DEBUG: BedGraph file size: $(stat -c%s "$BEDGRAPH_FILE") bytes"
        echo "DEBUG: Number of intervals: $(wc -l < "$BEDGRAPH_FILE")"
    else
        echo "ERROR: Failed to create bedGraph"
        exit 1
    fi
    
    # Convert bedGraph to BigWig
    echo "DEBUG: Converting bedGraph to BigWig..."
    bedGraphToBigWig "$BEDGRAPH_FILE" "$CHROM_SIZES" "$BIGWIG_FILE"
    
    if [[ $? -eq 0 ]]; then
        echo "DEBUG: BigWig created successfully"
        echo "DEBUG: BigWig file size: $(stat -c%s "$BIGWIG_FILE") bytes"
        
        # Clean up intermediate bedGraph file
        rm -f "$BEDGRAPH_FILE"
        echo "DEBUG: Cleaned up intermediate bedGraph file"
    else
        echo "ERROR: Failed to convert bedGraph to BigWig"
        exit 1
    fi
else
    echo "DEBUG: Skipping BigWig creation, file already exists: $BIGWIG_FILE"
fi

# Create normalized BigWig (optional - RPM normalized)
NORMALIZED_BIGWIG="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage_rpm.bw"
NORMALIZED_BEDGRAPH="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage_rpm.bedGraph"

if [[ ! -f "$NORMALIZED_BIGWIG" ]]; then
    echo "DEBUG: Creating RPM-normalized BigWig..."
    
    # Count total reads for normalization
    TOTAL_READS=$(zcat "$READS_FILE" | wc -l)
    SCALE_FACTOR=$(echo "scale=6; 1000000 / $TOTAL_READS" | bc)
    
    echo "DEBUG: Total reads: $(printf "%'d" $TOTAL_READS)"
    echo "DEBUG: RPM scale factor: $SCALE_FACTOR"
    
    # Create RPM-normalized bedGraph
    zcat "$READS_FILE" | \
        bedtools genomecov -i stdin -g "$CHROM_SIZES" -bg -scale $SCALE_FACTOR > "$NORMALIZED_BEDGRAPH"
    
    if [[ $? -eq 0 ]]; then
        echo "DEBUG: RPM-normalized bedGraph created"
        
        # Convert to BigWig
        bedGraphToBigWig "$NORMALIZED_BEDGRAPH" "$CHROM_SIZES" "$NORMALIZED_BIGWIG"
        
        if [[ $? -eq 0 ]]; then
            echo "DEBUG: RPM-normalized BigWig created successfully"
            rm -f "$NORMALIZED_BEDGRAPH"
            echo "DEBUG: Cleaned up intermediate RPM bedGraph file"
        else
            echo "ERROR: Failed to convert RPM bedGraph to BigWig"
            exit 1
        fi
    else
        echo "ERROR: Failed to create RPM-normalized bedGraph"
        exit 1
    fi
else
    echo "DEBUG: Skipping RPM-normalized BigWig creation, file already exists: $NORMALIZED_BIGWIG"
fi

# Create summary statistics
echo "DEBUG: Creating BigWig summary statistics..."
cat > "$OUTPUT_DIR/bigwig/${SAMPLE}_bigwig_info.txt" << EOF
SAMPLE=$SAMPLE
BIGWIG_CREATED_AT=$(date)
INPUT_READS_FILE=$READS_FILE
TOTAL_READS=$TOTAL_READS
RPM_SCALE_FACTOR=$SCALE_FACTOR
RAW_BIGWIG_FILE=${SAMPLE}_coverage.bw
NORMALIZED_BIGWIG_FILE=${SAMPLE}_coverage_rpm.bw
RAW_BIGWIG_SIZE=$(stat -c%s "$BIGWIG_FILE")
NORMALIZED_BIGWIG_SIZE=$(stat -c%s "$NORMALIZED_BIGWIG")
EOF

echo "BigWig files created:"
echo "  - Raw coverage: $OUTPUT_DIR/bigwig/${SAMPLE}_coverage.bw"
echo "  - RPM normalized: $OUTPUT_DIR/bigwig/${SAMPLE}_coverage_rpm.bw"
echo "  - Summary info: $OUTPUT_DIR/bigwig/${SAMPLE}_bigwig_info.txt"

echo "========================================="
echo "Step 7 complete for $SAMPLE"
echo "Created BigWig files for $(printf "%'d" $TOTAL_READS) reads"
echo "End time: $(date)"
echo "========================================="