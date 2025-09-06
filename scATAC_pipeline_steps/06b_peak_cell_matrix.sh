#!/bin/bash
#SBATCH --job-name=peak_cell_matrix_v2
#SBATCH --output=logs/06b_peak_cell_matrix_%a.out
#SBATCH --error=logs/06b_peak_cell_matrix_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# scATAC-seq Peak-Cell Matrix Generation (Optimized Version)
# Author: Pipeline Troubleshooting Analysis
# Date: September 2025
# Description: Generates peak-by-cell count matrices in 10X Genomics format
#              Optimized for large datasets with robust error handling

set -euo pipefail

# ==========================================
# CONFIGURATION & SETUP
# ==========================================

# Sample configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# Directory paths
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data"
OUTPUT_DIR="$BASE_DIR/chromap_final_output"
SCRIPT_DIR="$BASE_DIR/scATAC_pipeline_steps"
TEMP_DIR="/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# Create temporary directory
mkdir -p "$TEMP_DIR"
trap "rm -rf $TEMP_DIR" EXIT

# Set up conda environment
echo "Setting up conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

# ==========================================
# LOGGING & PROGRESS FUNCTIONS
# ==========================================

log_info() {
    echo "[INFO $(date '+%H:%M:%S')] $1" | tee -a "$SCRIPT_DIR/logs/06b_${SAMPLE}_progress.log"
}

log_error() {
    echo "[ERROR $(date '+%H:%M:%S')] $1" | tee -a "$SCRIPT_DIR/logs/06b_${SAMPLE}_progress.log"
}

log_progress() {
    local stage=$1 percent=$2
    echo "[PROGRESS $(date '+%H:%M:%S')] [$stage] $percent% complete" | tee -a "$SCRIPT_DIR/logs/06b_${SAMPLE}_progress.log"
}

report_resource_usage() {
    if command -v free &> /dev/null; then
        echo "Memory usage: $(free -h | grep '^Mem:' | awk '{print $3"/"$2}')"
    fi
    echo "Disk usage (temp): $(du -sh $TEMP_DIR 2>/dev/null || echo 'N/A')"
}

# ==========================================
# INPUT VALIDATION FUNCTIONS
# ==========================================

validate_file() {
    local file=$1 min_size=${2:-1}
    local description=${3:-"File"}
    
    if [[ ! -f "$file" ]]; then
        log_error "$description not found: $file"
        return 1
    fi
    
    local size=$(stat -c%s "$file" 2>/dev/null || echo 0)
    if [[ $size -lt $min_size ]]; then
        log_error "$description is too small ($size bytes): $file"
        return 1
    fi
    
    log_info "$description validated: $file ($size bytes)"
    return 0
}

validate_prerequisites() {
    log_info "Validating prerequisites for sample: $SAMPLE"
    
    # Define required files
    local peaks_file="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"
    local reads_file="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
    local genome_sizes="$OUTPUT_DIR/qc/genome.chrom.sizes"
    
    # Validate each required file
    validate_file "$peaks_file" 1000 "Sorted peaks file" || exit 1
    validate_file "$reads_file" 1000000 "Reads file" || exit 1
    validate_file "$genome_sizes" 100 "Genome sizes file" || exit 1
    
    # Check if previous matrix exists and warn
    local matrix_dir="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
    if [[ -d "$matrix_dir" ]]; then
        log_info "Previous matrix directory exists, will be overwritten: $matrix_dir"
        rm -rf "$matrix_dir"
    fi
    
    # Validate tools
    if ! command -v bedtools &> /dev/null; then
        log_error "bedtools not found. Please install bedtools."
        exit 1
    fi
    
    if ! command -v sort &> /dev/null; then
        log_error "sort command not found."
        exit 1
    fi
    
    log_info "All prerequisites validated successfully"
    return 0
}

# ==========================================
# CORE PROCESSING FUNCTIONS
# ==========================================

estimate_resources() {
    local peaks_file="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"
    local reads_file="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
    
    local peak_count=$(wc -l < "$peaks_file")
    local reads_size=$(stat -c%s "$reads_file")
    
    log_info "Dataset estimation:"
    log_info "  Peak count: $(printf "%'d" $peak_count)"
    log_info "  Reads file size: $(numfmt --to=iec-i --suffix=B $reads_size)"
    log_info "  Estimated overlap records: $(printf "%'d" $((peak_count * reads_size / 100000000)))"
}

run_bedtools_intersect() {
    local peaks_file="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"
    local reads_file="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
    local genome_sizes="$OUTPUT_DIR/qc/genome.chrom.sizes"
    local output_file="$TEMP_DIR/${SAMPLE}_intersect_raw.txt"
    
    log_info "Starting bedtools intersect analysis..."
    log_progress "bedtools_intersect" 0
    report_resource_usage
    
    # Run bedtools intersect with progress monitoring
    {
        zcat "$reads_file" | \
        sort -k1,1 -k2,2n -S32G --parallel=$SLURM_CPUS_PER_TASK --temporary-directory="$TEMP_DIR" | \
        bedtools intersect \
            -a "$peaks_file" \
            -b stdin \
            -wo -sorted -g "$genome_sizes"
    } > "$output_file" &
    
    local intersect_pid=$!
    
    # Monitor progress
    while kill -0 $intersect_pid 2>/dev/null; do
        sleep 30
        if [[ -f "$output_file" ]]; then
            local current_lines=$(wc -l < "$output_file" 2>/dev/null || echo 0)
            log_info "Bedtools progress: $(printf "%'d" $current_lines) overlaps found"
        fi
        report_resource_usage
    done
    
    wait $intersect_pid
    local exit_code=$?
    
    if [[ $exit_code -ne 0 ]]; then
        log_error "Bedtools intersect failed with exit code: $exit_code"
        return 1
    fi
    
    if [[ ! -s "$output_file" ]]; then
        log_error "Bedtools intersect produced empty output"
        return 1
    fi
    
    local overlap_count=$(wc -l < "$output_file")
    log_info "Bedtools intersect completed: $(printf "%'d" $overlap_count) overlaps found"
    log_progress "bedtools_intersect" 100
    
    echo "$output_file"
    return 0
}

process_overlaps() {
    local intersect_file=$1
    local output_file="$TEMP_DIR/${SAMPLE}_peak_barcode_counts.tsv.gz"
    
    log_info "Processing overlap data to generate peak-barcode counts..."
    log_progress "process_overlaps" 0
    report_resource_usage
    
    # Extract peak ID and barcode, count occurrences
    {
        cat "$intersect_file" | \
        awk -v OFS='\t' '{print $4, $8}' | \
        sort -k1,1 -k2,2 -S64G --parallel=$SLURM_CPUS_PER_TASK --temporary-directory="$TEMP_DIR" | \
        uniq -c | \
        awk -v OFS='\t' '{print $2, $3, $1}' | \
        gzip -c
    } > "$output_file" &
    
    local process_pid=$!
    
    # Monitor progress
    while kill -0 $process_pid 2>/dev/null; do
        sleep 15
        report_resource_usage
    done
    
    wait $process_pid
    local exit_code=$?
    
    if [[ $exit_code -ne 0 ]]; then
        log_error "Overlap processing failed with exit code: $exit_code"
        return 1
    fi
    
    if [[ ! -s "$output_file" ]]; then
        log_error "Overlap processing produced empty output"
        return 1
    fi
    
    local entry_count=$(zcat "$output_file" | wc -l)
    log_info "Overlap processing completed: $(printf "%'d" $entry_count) peak-barcode pairs"
    log_progress "process_overlaps" 100
    
    echo "$output_file"
    return 0
}

create_10x_matrix() {
    local peak_barcode_file=$1
    local peaks_file="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"
    local matrix_dir="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
    
    log_info "Creating 10X Genomics format matrix..."
    log_progress "create_matrix" 0
    
    # Create output directory
    mkdir -p "$matrix_dir"
    
    # 1. Create features.tsv (peaks)
    log_info "Generating features.tsv..."
    awk -v OFS='\t' '{print $4, $4, "Peaks"}' "$peaks_file" > "$matrix_dir/features.tsv"
    local feature_count=$(wc -l < "$matrix_dir/features.tsv")
    log_info "Created features file with $(printf "%'d" $feature_count) peaks"
    
    # 2. Create barcodes.tsv
    log_info "Generating barcodes.tsv..."
    zcat "$peak_barcode_file" | cut -f2 | sort -u > "$matrix_dir/barcodes.tsv"
    local barcode_count=$(wc -l < "$matrix_dir/barcodes.tsv")
    log_info "Created barcodes file with $(printf "%'d" $barcode_count) cells"
    
    if [[ $barcode_count -eq 0 ]]; then
        log_error "No valid barcodes found"
        return 1
    fi
    
    log_progress "create_matrix" 50
    
    # 3. Create matrix.mtx with indexed coordinates
    log_info "Generating matrix.mtx..."
    local temp_matrix="$TEMP_DIR/matrix.mtx.tmp"
    
    # Create AWK script for efficient matrix generation
    awk -v peak_file="$matrix_dir/features.tsv" -v barcode_file="$matrix_dir/barcodes.tsv" '
    BEGIN {
        OFS="\t"
        # Read peak IDs into array
        while ((getline line < peak_file) > 0) {
            split(line, fields, "\t")
            peak_map[fields[1]] = NR
        }
        close(peak_file)
        
        # Read barcode IDs into array
        FNR_barcode = 0
        while ((getline line < barcode_file) > 0) {
            FNR_barcode++
            barcode_map[line] = FNR_barcode
        }
        close(barcode_file)
    }
    {
        peak_id = $1
        barcode_id = $2
        count = $3
        
        peak_idx = peak_map[peak_id]
        barcode_idx = barcode_map[barcode_id]
        
        if (peak_idx > 0 && barcode_idx > 0) {
            print peak_idx, barcode_idx, count
        }
    }' <(zcat "$peak_barcode_file") > "$temp_matrix"
    
    local entry_count=$(wc -l < "$temp_matrix")
    log_info "Generated matrix with $(printf "%'d" $entry_count) entries"
    
    # Create final MTX format with header
    {
        echo "%%MatrixMarket matrix coordinate integer general"
        echo "% Generated by scATAC pipeline step 6b on $(date)"
        echo "$feature_count $barcode_count $entry_count"
        cat "$temp_matrix"
    } > "$matrix_dir/matrix.mtx"
    
    log_progress "create_matrix" 75
    
    # 4. Compress files
    log_info "Compressing matrix files..."
    gzip -f "$matrix_dir/matrix.mtx"
    gzip -f "$matrix_dir/features.tsv"
    gzip -f "$matrix_dir/barcodes.tsv"
    
    # 5. Calculate and save statistics
    local total_possible_entries=$((feature_count * barcode_count))
    local sparsity_percent=0
    if [[ $total_possible_entries -gt 0 ]]; then
        sparsity_percent=$(echo "scale=2; (1 - $entry_count / $total_possible_entries) * 100" | bc -l 2>/dev/null || echo "N/A")
    fi
    
    # Create comprehensive matrix info
    cat > "$matrix_dir/matrix_info.txt" << EOF
SAMPLE=$SAMPLE
TOTAL_PEAKS=$feature_count
TOTAL_BARCODES=$barcode_count
TOTAL_ENTRIES=$entry_count
SPARSITY=${sparsity_percent}%
MATRIX_CREATED_AT=$(date)
PIPELINE_VERSION=6b_optimized
SLURM_JOB_ID=$SLURM_JOB_ID
RESOURCES_USED=CPUs:$SLURM_CPUS_PER_TASK,MEM:${SLURM_MEM_PER_NODE}MB
INPUT_PEAKS_FILE=$peaks_file
INPUT_READS_FILE=$OUTPUT_DIR/${SAMPLE}_reads.bed.gz
PROCESSING_TIME=$(date -d@$SECONDS -u +%H:%M:%S)
EOF
    
    # Copy peak-barcode counts for downstream analysis
    cp "$peak_barcode_file" "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz"
    
    log_info "Matrix statistics:"
    log_info "  Peaks: $(printf "%'d" $feature_count)"
    log_info "  Cells: $(printf "%'d" $barcode_count)"  
    log_info "  Entries: $(printf "%'d" $entry_count)"
    log_info "  Sparsity: $sparsity_percent%"
    
    log_progress "create_matrix" 100
    return 0
}

validate_output() {
    local matrix_dir="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
    
    log_info "Validating output files..."
    
    # Check required files exist and are non-empty
    local required_files=(
        "$matrix_dir/matrix.mtx.gz"
        "$matrix_dir/features.tsv.gz" 
        "$matrix_dir/barcodes.tsv.gz"
        "$matrix_dir/matrix_info.txt"
    )
    
    for file in "${required_files[@]}"; do
        if [[ ! -f "$file" ]] || [[ ! -s "$file" ]]; then
            log_error "Output validation failed: $file"
            return 1
        fi
    done
    
    # Validate matrix dimensions consistency
    local info_peaks=$(grep "TOTAL_PEAKS" "$matrix_dir/matrix_info.txt" | cut -d'=' -f2)
    local info_barcodes=$(grep "TOTAL_BARCODES" "$matrix_dir/matrix_info.txt" | cut -d'=' -f2)
    local info_entries=$(grep "TOTAL_ENTRIES" "$matrix_dir/matrix_info.txt" | cut -d'=' -f2)
    
    if [[ $info_entries -eq 0 ]]; then
        log_error "Matrix has zero entries - this indicates a processing failure"
        return 1
    fi
    
    log_info "Output validation successful:"
    log_info "  Matrix dimensions: $info_peaks × $info_barcodes"
    log_info "  Non-zero entries: $(printf "%'d" $info_entries)"
    
    return 0
}

# ==========================================
# MAIN EXECUTION
# ==========================================

main() {
    local start_time=$(date +%s)
    
    echo "========================================="
    echo "Step 6b: Optimized Peak-Cell Matrix Generation"
    echo "Sample: $SAMPLE"
    echo "Start time: $(date)"
    echo "SLURM Job ID: $SLURM_JOB_ID"
    echo "Resources: ${SLURM_CPUS_PER_TASK} CPUs, ${SLURM_MEM_PER_NODE}MB RAM"
    echo "========================================="
    
    # Step 1: Validate prerequisites
    log_progress "validation" 0
    validate_prerequisites || exit 1
    log_progress "validation" 100
    
    # Step 2: Estimate resource requirements
    estimate_resources
    
    # Step 3: Run bedtools intersect
    log_progress "overall" 10
    local intersect_file
    intersect_file=$(run_bedtools_intersect) || exit 1
    log_progress "overall" 40
    
    # Step 4: Process overlaps to create peak-barcode counts
    local peak_barcode_file
    peak_barcode_file=$(process_overlaps "$intersect_file") || exit 1
    log_progress "overall" 70
    
    # Step 5: Create 10X format matrix
    create_10x_matrix "$peak_barcode_file" || exit 1
    log_progress "overall" 90
    
    # Step 6: Validate output
    validate_output || exit 1
    log_progress "overall" 100
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    echo "========================================="
    echo "Step 6b completed successfully for $SAMPLE"
    echo "Processing time: $(date -d@$duration -u +%H:%M:%S)"
    echo "Output directory: $OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
    echo "Files generated:"
    echo "  - matrix.mtx.gz (sparse matrix)"
    echo "  - features.tsv.gz (peak annotations)"
    echo "  - barcodes.tsv.gz (cell barcodes)"
    echo "  - matrix_info.txt (metadata)"
    echo "End time: $(date)"
    echo "========================================="
}

# Execute main function
main "$@"