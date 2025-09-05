#!/bin/bash
#SBATCH --job-name=qc_metrics
#SBATCH --output=logs/08_qc_metrics_%a.out
#SBATCH --error=logs/08_qc_metrics_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 8: QC Metrics Analysis for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create QC output directory
mkdir -p "$OUTPUT_DIR/qc_metrics"

# Check prerequisites
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"
CHROM_SIZES="$OUTPUT_DIR/mm10.chrom.sizes"

if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    exit 1
fi

if [[ ! -f "$PEAKS_FILE" ]]; then
    echo "ERROR: Peaks file not found: $PEAKS_FILE"
    echo "Please run step 5 (05_call_peaks.sh) first"
    exit 1
fi

# Check for required tools
echo "DEBUG: Checking for required tools..."
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found"
    exit 1
fi

echo "DEBUG: Tools verified successfully"

# 1. Fragment length distribution
echo "DEBUG: Calculating fragment length distribution..."
zcat "$FRAGMENTS_FILE" | \
    awk '{print $3-$2}' | \
    sort -n | \
    uniq -c | \
    awk '{print $2"\t"$1}' > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragment_lengths.txt"

# Calculate fragment length statistics
FRAGMENT_STATS="$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragment_stats.txt"
zcat "$FRAGMENTS_FILE" | \
    awk '{len=$3-$2; sum+=len; sumsq+=len*len; if(NR==1){min=max=len} if(len<min){min=len} if(len>max){max=len}} 
         END{mean=sum/NR; variance=(sumsq-sum*sum/NR)/(NR-1); sd=sqrt(variance); 
         print "Total_fragments\t"NR"\nMean_length\t"mean"\nMedian_length\t"median"\nSD_length\t"sd"\nMin_length\t"min"\nMax_length\t"max}' > "$FRAGMENT_STATS"

echo "DEBUG: Fragment length analysis completed"

# 2. TSS enrichment calculation
echo "DEBUG: Calculating TSS enrichment..."
# Download TSS regions if not available (mouse mm10)
TSS_FILE="$OUTPUT_DIR/qc_metrics/mm10_tss.bed"
if [[ ! -f "$TSS_FILE" ]]; then
    echo "DEBUG: Downloading comprehensive mm10 TSS file from UCSC..."
    # Download and process RefSeq gene annotations for mm10 from UCSC
    wget -qO- "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz" | \
        gunzip -c | \
        awk 'BEGIN{OFS="\t"} {if($4=="+"){print $3, $5-1, $5, $13, "0", $4} else {print $3, $6-1, $6, $13, "0", $4}}' | \
        sort -k1,1 -k2,2n > "$TSS_FILE"
    echo "DEBUG: Comprehensive TSS file created successfully."
fi

# Calculate TSS enrichment (fragments overlapping TSS +/- 2kb)
bedtools slop -i "$TSS_FILE" -g "$CHROM_SIZES" -b 2000 | \
    bedtools intersect -a "$FRAGMENTS_FILE" -b stdin -c > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_overlap.txt"

# Calculate enrichment score
TSS_ENRICHMENT=$(awk 'BEGIN{tss=0; total=0} {total++; if($5>0) tss++} END{print tss/total*100}' "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_overlap.txt")
echo "TSS_enrichment_percent\t$TSS_ENRICHMENT" >> "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_stats.txt"

echo "DEBUG: TSS enrichment calculated: $TSS_ENRICHMENT%"

# 3. Peak-in-promoter analysis
echo "DEBUG: Analyzing peaks in promoter regions..."
if [[ -f "$TSS_FILE" ]]; then
    # Create promoter regions (TSS +/- 1kb)
    bedtools slop -i "$TSS_FILE" -g "$CHROM_SIZES" -b 1000 > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoters.bed"
    
    # Count peaks in promoters
    PEAKS_IN_PROMOTERS=$(bedtools intersect -a "$PEAKS_FILE" -b "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoters.bed" -u | wc -l)
    TOTAL_PEAKS=$(wc -l < "$PEAKS_FILE")
    PROMOTER_PERCENT=$(awk -v p="$PEAKS_IN_PROMOTERS" -v t="$TOTAL_PEAKS" 'BEGIN {printf "%.2f", p/t*100}')
    
    echo "Peaks_in_promoters\t$PEAKS_IN_PROMOTERS" > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt"
    echo "Total_peaks\t$TOTAL_PEAKS" >> "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt"
    echo "Promoter_percentage\t$PROMOTER_PERCENT" >> "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt"
    
    echo "DEBUG: Promoter analysis completed: $PROMOTER_PERCENT% peaks in promoters"
fi

# 4. Cell quality metrics
echo "DEBUG: Calculating per-cell quality metrics..."
# Calculate fragments per cell
zcat "$FRAGMENTS_FILE" | \
    cut -f4 | \
    sort | \
    uniq -c | \
    awk '{print $2"\t"$1}' > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragments_per_cell.txt"

# Calculate peaks per cell (if peak-cell matrix exists)
PEAK_BC_MATRIX="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
if [[ -d "$PEAK_BC_MATRIX" ]]; then
    if [[ -f "$PEAK_BC_MATRIX/matrix.mtx.gz" ]]; then
        echo "DEBUG: Calculating peaks per cell from matrix..."
        # Extract non-zero entries per barcode
        zcat "$PEAK_BC_MATRIX/matrix.mtx.gz" | \
            tail -n +4 | \
            awk '{count[$2]++} END{for(bc in count) print bc"\t"count[bc]}' | \
            sort -n > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_peaks_per_cell_idx.txt"
    fi
fi

# 5. Library complexity
echo "DEBUG: Calculating library complexity..."
TOTAL_FRAGMENTS=$(zcat "$FRAGMENTS_FILE" | wc -l)
UNIQUE_FRAGMENTS=$(zcat "$FRAGMENTS_FILE" | cut -f1-3 | sort -u | wc -l)
COMPLEXITY=$(echo "scale=4; $UNIQUE_FRAGMENTS / $TOTAL_FRAGMENTS" | bc)

cat > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_library_complexity.txt" << EOF
Total_fragments	$TOTAL_FRAGMENTS
Unique_fragments	$UNIQUE_FRAGMENTS
Library_complexity	$COMPLEXITY
EOF

echo "DEBUG: Library complexity: $COMPLEXITY"

# 6. Create comprehensive QC summary
echo "DEBUG: Creating comprehensive QC summary..."
cat > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_qc_summary.txt" << EOF
QC_METRICS_SUMMARY
Sample: $SAMPLE
Analysis_date: $(date)

FRAGMENT_METRICS:
$(cat "$FRAGMENT_STATS")

TSS_METRICS:
$(cat "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_stats.txt" 2>/dev/null || echo "TSS_enrichment_percent	N/A")

PROMOTER_METRICS:
$(cat "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt" 2>/dev/null || echo "Promoter_analysis	Not_available")

LIBRARY_COMPLEXITY:
$(cat "$OUTPUT_DIR/qc_metrics/${SAMPLE}_library_complexity.txt")

CELL_COUNT:
Total_barcodes	$(wc -l < "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragments_per_cell.txt")

FILES_GENERATED:
- Fragment length distribution: ${SAMPLE}_fragment_lengths.txt
- Fragment statistics: ${SAMPLE}_fragment_stats.txt
- TSS enrichment: ${SAMPLE}_tss_stats.txt
- Promoter analysis: ${SAMPLE}_promoter_stats.txt
- Fragments per cell: ${SAMPLE}_fragments_per_cell.txt
- Library complexity: ${SAMPLE}_library_complexity.txt
- QC summary: ${SAMPLE}_qc_summary.txt
EOF

# 7. Generate R script for QC plotting (optional)
cat > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_plot_qc.R" << 'EOF'
#!/usr/bin/env Rscript
# QC Plotting Script for scATAC-seq data
# Usage: Rscript plot_qc.R sample_name

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  stop("Please provide sample name as argument")
}

sample_name <- args[1]
library(ggplot2)
library(dplyr)

# Set working directory to qc_metrics folder
setwd("qc_metrics")

# 1. Fragment length distribution
frag_lengths <- read.table(paste0(sample_name, "_fragment_lengths.txt"), 
                          col.names = c("Length", "Count"))

p1 <- ggplot(frag_lengths, aes(x = Length, y = Count)) +
  geom_line() +
  xlim(0, 1000) +
  labs(title = paste("Fragment Length Distribution -", sample_name),
       x = "Fragment Length (bp)", y = "Count") +
  theme_minimal()

ggsave(paste0(sample_name, "_fragment_length_dist.png"), p1, width = 8, height = 6)

# 2. Fragments per cell distribution
frags_per_cell <- read.table(paste0(sample_name, "_fragments_per_cell.txt"), 
                            col.names = c("Barcode", "Fragment_count"))

p2 <- ggplot(frags_per_cell, aes(x = Fragment_count)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_x_log10() +
  labs(title = paste("Fragments per Cell Distribution -", sample_name),
       x = "Fragments per Cell (log10)", y = "Number of Cells") +
  theme_minimal()

ggsave(paste0(sample_name, "_fragments_per_cell.png"), p2, width = 8, height = 6)

cat("QC plots saved successfully\n")
EOF

chmod +x "$OUTPUT_DIR/qc_metrics/${SAMPLE}_plot_qc.R"

echo "QC Analysis Summary:"
echo "===================="
echo "Total fragments: $(printf "%'d" $TOTAL_FRAGMENTS)"
echo "Unique fragments: $(printf "%'d" $UNIQUE_FRAGMENTS)"
echo "Library complexity: $COMPLEXITY"
echo "TSS enrichment: $TSS_ENRICHMENT%"
if [[ -f "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt" ]]; then
    echo "Promoter peaks: $PROMOTER_PERCENT%"
fi
echo "Total barcodes: $(wc -l < "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragments_per_cell.txt")"

echo "Output files created in: $OUTPUT_DIR/qc_metrics/"
echo "  - Comprehensive QC summary: ${SAMPLE}_qc_summary.txt"
echo "  - Fragment statistics: ${SAMPLE}_fragment_stats.txt"
echo "  - Per-cell metrics: ${SAMPLE}_fragments_per_cell.txt"
echo "  - R plotting script: ${SAMPLE}_plot_qc.R"

echo "========================================="
echo "Step 8 complete for $SAMPLE"
echo "QC metrics analysis finished"
echo "End time: $(date)"
echo "========================================="