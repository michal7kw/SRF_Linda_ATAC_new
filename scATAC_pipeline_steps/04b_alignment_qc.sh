#!/bin/bash
#SBATCH --job-name=alignment_qc
#SBATCH --output=logs/04b_alignment_qc_%a.out
#SBATCH --error=logs/04b_alignment_qc_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
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
echo "Step 4b: Alignment Quality Control for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create alignment QC output directory
mkdir -p "$OUTPUT_DIR/alignment_qc"

# Check prerequisites
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
CHROMAP_LOG="$OUTPUT_DIR/logs/${SAMPLE}_chromap.log"

if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    echo "Please run step 4 (04_chromap_alignment.sh) first"
    exit 1
fi

if [[ ! -f "$READS_FILE" ]]; then
    echo "ERROR: Reads file not found: $READS_FILE"
    echo "Please run step 4 (04_chromap_alignment.sh) first"
    exit 1
fi

echo "DEBUG: Prerequisites verified"

# Check for required tools
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found"
    exit 1
fi

echo "DEBUG: Required tools verified"

# 1. Parse chromap alignment statistics
echo "DEBUG: Parsing chromap alignment statistics..."
if [[ -f "$CHROMAP_LOG" ]]; then
    # Extract key alignment metrics from chromap log
    cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt" << EOF
# Chromap Alignment Statistics for $SAMPLE
# Generated: $(date)

EOF
    
    # Extract total reads
    TOTAL_READS=$(grep -i "total.*reads" "$CHROMAP_LOG" | head -1 | grep -oE '[0-9,]+' | tr -d ',' || echo "N/A")
    echo "Total_input_reads\t$TOTAL_READS" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    
    # Extract mapped reads
    MAPPED_READS=$(grep -i "mapped.*reads\|aligned.*reads" "$CHROMAP_LOG" | head -1 | grep -oE '[0-9,]+' | tr -d ',' || echo "N/A")
    echo "Mapped_reads\t$MAPPED_READS" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    
    # Calculate mapping rate
    if [[ "$TOTAL_READS" != "N/A" && "$MAPPED_READS" != "N/A" ]]; then
        MAPPING_RATE=$(echo "scale=2; $MAPPED_READS * 100 / $TOTAL_READS" | bc)
        echo "Mapping_rate_percent\t$MAPPING_RATE" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    else
        echo "Mapping_rate_percent\tN/A" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    fi
    
    # Extract duplicate rate if available
    DUPLICATE_RATE=$(grep -i "duplicate" "$CHROMAP_LOG" | grep -oE '[0-9.]+%' | head -1 || echo "N/A")
    echo "Duplicate_rate\t$DUPLICATE_RATE" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    
    echo "DEBUG: Chromap statistics extracted"
else
    echo "WARNING: Chromap log file not found: $CHROMAP_LOG"
fi

# 2. Fragment statistics analysis
echo "DEBUG: Analyzing fragment statistics..."

# Basic fragment counts
TOTAL_FRAGMENTS=$(zcat "$FRAGMENTS_FILE" | wc -l)
UNIQUE_BARCODES=$(zcat "$FRAGMENTS_FILE" | cut -f4 | sort -u | wc -l)

echo "Total_fragments\t$TOTAL_FRAGMENTS" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_summary.txt"
echo "Unique_barcodes\t$UNIQUE_BARCODES" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_summary.txt"

# Calculate fragments per barcode statistics
echo "DEBUG: Calculating per-barcode fragment statistics..."
zcat "$FRAGMENTS_FILE" | \
    cut -f4 | \
    sort | \
    uniq -c | \
    awk '{print $1}' | \
    sort -n > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt"

# Calculate summary statistics
FRAGMENT_STATS="$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_distribution.txt"
awk 'BEGIN{sum=0; count=0; min=""; max=""} 
     {
         sum+=$1; count++; 
         if(min=="" || $1<min) min=$1; 
         if(max=="" || $1>max) max=$1;
         values[count]=$1
     } 
     END{
         mean=sum/count
         asort(values)
         median=(count%2==1) ? values[int(count/2)+1] : (values[int(count/2)] + values[int(count/2)+1])/2
         print "Mean_fragments_per_barcode\t"mean
         print "Median_fragments_per_barcode\t"median
         print "Min_fragments_per_barcode\t"min
         print "Max_fragments_per_barcode\t"max
         print "Total_barcodes_with_fragments\t"count
     }' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" > "$FRAGMENT_STATS"

echo "DEBUG: Fragment distribution analysis completed"

# 3. Fragment length distribution analysis
echo "DEBUG: Analyzing fragment length distribution..."
zcat "$FRAGMENTS_FILE" | \
    awk '{print $3-$2}' | \
    sort -n | \
    uniq -c | \
    awk '{print $2"\t"$1}' > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_lengths_detailed.txt"

# Calculate fragment length statistics
FRAGMENT_LENGTH_STATS="$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_length_stats.txt"
zcat "$FRAGMENTS_FILE" | \
    awk '{len=$3-$2; sum+=len; sumsq+=len*len; lengths[NR]=len} 
         END{
             mean=sum/NR; 
             for(i=1;i<=NR;i++) variance+=(lengths[i]-mean)^2
             variance=variance/(NR-1)
             sd=sqrt(variance)
             asort(lengths)
             median=(NR%2==1) ? lengths[int(NR/2)+1] : (lengths[int(NR/2)] + lengths[int(NR/2)+1])/2
             q25=lengths[int(NR*0.25)]
             q75=lengths[int(NR*0.75)]
             
             print "Mean_fragment_length\t"mean
             print "Median_fragment_length\t"median
             print "SD_fragment_length\t"sd
             print "Q25_fragment_length\t"q25
             print "Q75_fragment_length\t"q75
             
             # Count fragments in key size ranges
             mono_nucleosome=0; di_nucleosome=0; long_fragments=0
             for(i=1;i<=NR;i++) {
                 if(lengths[i]<=150) mono_nucleosome++
                 else if(lengths[i]<=300) di_nucleosome++
                 else long_fragments++
             }
             
             print "Mono_nucleosome_fragments_<=150bp\t"mono_nucleosome
             print "Di_nucleosome_fragments_150-300bp\t"di_nucleosome
             print "Long_fragments_>300bp\t"long_fragments
             print "Mono_nucleosome_percentage\t"(mono_nucleosome*100/NR)
             print "Di_nucleosome_percentage\t"(di_nucleosome*100/NR)
         }' > "$FRAGMENT_LENGTH_STATS"

echo "DEBUG: Fragment length analysis completed"

# 4. Chromosome distribution analysis
echo "DEBUG: Analyzing chromosome distribution..."
zcat "$FRAGMENTS_FILE" | \
    cut -f1 | \
    sort | \
    uniq -c | \
    awk '{print $2"\t"$1}' | \
    sort -k1,1V > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_chromosome_distribution.txt"

# Calculate mitochondrial fragment percentage
MITO_FRAGMENTS=$(grep -E "^chrM|^MT" "$OUTPUT_DIR/alignment_qc/${SAMPLE}_chromosome_distribution.txt" | awk '{sum+=$2} END{print sum+0}')
MITO_PERCENTAGE=$(echo "scale=2; $MITO_FRAGMENTS * 100 / $TOTAL_FRAGMENTS" | bc)

echo "Mitochondrial_fragments\t$MITO_FRAGMENTS" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt"
echo "Mitochondrial_percentage\t$MITO_PERCENTAGE" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt"

echo "DEBUG: Mitochondrial percentage: $MITO_PERCENTAGE%"

# 5. Read quality assessment
echo "DEBUG: Analyzing read quality from BED file..."
if [[ -f "$READS_FILE" ]]; then
    TOTAL_READS_BED=$(zcat "$READS_FILE" | wc -l)
    UNIQUE_READ_POSITIONS=$(zcat "$READS_FILE" | cut -f1-3 | sort -u | wc -l)
    
    echo "Total_reads_in_bed\t$TOTAL_READS_BED" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt"
    echo "Unique_read_positions\t$UNIQUE_READ_POSITIONS" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt"
    
    # Calculate PCR duplicate rate
    if [[ $TOTAL_READS_BED -gt 0 ]]; then
        DUPLICATE_RATE_CALC=$(echo "scale=4; (1 - $UNIQUE_READ_POSITIONS / $TOTAL_READS_BED) * 100" | bc)
        echo "Estimated_PCR_duplicate_rate_percent\t$DUPLICATE_RATE_CALC" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt"
    fi
fi

# 6. Barcode quality assessment
echo "DEBUG: Assessing barcode quality..."
# Barcode complexity
BARCODE_COMPLEXITY=$(echo "scale=4; $UNIQUE_BARCODES / $TOTAL_FRAGMENTS" | bc)
echo "Barcode_complexity\t$BARCODE_COMPLEXITY" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"

# Distribution of fragments per barcode
BARCODES_GT_100=$(awk '$1>=100' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" | wc -l)
BARCODES_GT_500=$(awk '$1>=500' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" | wc -l)
BARCODES_GT_1000=$(awk '$1>=1000' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" | wc -l)

echo "Barcodes_with_100plus_fragments\t$BARCODES_GT_100" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"
echo "Barcodes_with_500plus_fragments\t$BARCODES_GT_500" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"
echo "Barcodes_with_1000plus_fragments\t$BARCODES_GT_1000" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"

# 7. Create comprehensive alignment QC summary
echo "DEBUG: Creating comprehensive QC summary..."
cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_qc_summary.txt" << EOF
# Alignment Quality Control Summary for $SAMPLE
# Generated: $(date)
# Analysis performed before peak calling

=== ALIGNMENT STATISTICS ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt" 2>/dev/null | grep -v "^#" || echo "Alignment_stats\tNot_available")

=== FRAGMENT SUMMARY ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_summary.txt")

=== FRAGMENT DISTRIBUTION ===
$(cat "$FRAGMENT_STATS")

=== FRAGMENT LENGTH STATISTICS ===
$(cat "$FRAGMENT_LENGTH_STATS")

=== MITOCHONDRIAL CONTENT ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt")

=== READ QUALITY ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt" 2>/dev/null || echo "Read_stats\tNot_available")

=== BARCODE QUALITY ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt")
EOF

# 8. Generate R script for alignment QC plots
echo "DEBUG: Creating R plotting script..."
cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_plot_alignment_qc.R" << 'EOF'
#!/usr/bin/env Rscript
# Alignment QC Plotting Script

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  stop("Please provide sample name as argument")
}

sample_name <- args[1]

suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
})

setwd("alignment_qc")

# 1. Fragment length distribution
if(file.exists(paste0(sample_name, "_fragment_lengths_detailed.txt"))) {
    frag_lengths <- read.table(paste0(sample_name, "_fragment_lengths_detailed.txt"), 
                              col.names = c("Length", "Count"))
    
    p1 <- ggplot(frag_lengths %>% filter(Length <= 1000), aes(x = Length, y = Count)) +
        geom_line(color = "blue", alpha = 0.8) +
        geom_vline(xintercept = 147, color = "red", linetype = "dashed", alpha = 0.7) +
        geom_vline(xintercept = 294, color = "orange", linetype = "dashed", alpha = 0.7) +
        annotate("text", x = 147, y = max(frag_lengths$Count[frag_lengths$Length <= 1000]) * 0.8, 
                 label = "Mono-nucleosome", angle = 90, vjust = -0.5, size = 3) +
        annotate("text", x = 294, y = max(frag_lengths$Count[frag_lengths$Length <= 1000]) * 0.8, 
                 label = "Di-nucleosome", angle = 90, vjust = -0.5, size = 3) +
        labs(title = paste("Fragment Length Distribution -", sample_name),
             x = "Fragment Length (bp)", y = "Count") +
        theme_minimal()
    
    ggsave(paste0(sample_name, "_fragment_length_distribution.png"), p1, width = 10, height = 6, dpi = 300)
}

# 2. Fragments per barcode distribution
if(file.exists(paste0(sample_name, "_fragments_per_barcode_raw.txt"))) {
    frags_per_bc <- read.table(paste0(sample_name, "_fragments_per_barcode_raw.txt"), 
                              col.names = "Fragment_count")
    
    p2 <- ggplot(frags_per_bc, aes(x = Fragment_count)) +
        geom_histogram(bins = 50, alpha = 0.7, fill = "skyblue") +
        scale_x_log10() +
        geom_vline(xintercept = 100, color = "red", linetype = "dashed") +
        geom_vline(xintercept = 500, color = "orange", linetype = "dashed") +
        geom_vline(xintercept = 1000, color = "green", linetype = "dashed") +
        labs(title = paste("Fragments per Barcode Distribution -", sample_name),
             x = "Fragments per Barcode (log10)", y = "Number of Barcodes") +
        theme_minimal()
    
    ggsave(paste0(sample_name, "_fragments_per_barcode.png"), p2, width = 10, height = 6, dpi = 300)
}

# 3. Chromosome distribution
if(file.exists(paste0(sample_name, "_chromosome_distribution.txt"))) {
    chrom_dist <- read.table(paste0(sample_name, "_chromosome_distribution.txt"), 
                            col.names = c("Chromosome", "Count"))
    
    # Filter for main chromosomes and sort
    main_chroms <- chrom_dist %>% 
        filter(grepl("^chr[0-9XY]+$", Chromosome)) %>%
        mutate(Chromosome = factor(Chromosome, 
                                  levels = paste0("chr", c(1:19, "X", "Y"))))
    
    p3 <- ggplot(main_chroms, aes(x = Chromosome, y = Count)) +
        geom_bar(stat = "identity", alpha = 0.7, fill = "lightcoral") +
        labs(title = paste("Fragment Distribution by Chromosome -", sample_name),
             x = "Chromosome", y = "Fragment Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0(sample_name, "_chromosome_distribution.png"), p3, width = 12, height = 6, dpi = 300)
}

cat("Alignment QC plots generated successfully!\n")
EOF

chmod +x "$OUTPUT_DIR/alignment_qc/${SAMPLE}_plot_alignment_qc.R"

# 9. Run R plotting if possible
if command -v Rscript &> /dev/null; then
    echo "DEBUG: Generating alignment QC plots..."
    cd "$OUTPUT_DIR"
    Rscript "alignment_qc/${SAMPLE}_plot_alignment_qc.R" "$SAMPLE" || echo "WARNING: R plotting failed"
    cd - > /dev/null
else
    echo "WARNING: Rscript not available. Plots not generated."
fi

# 10. Generate quality assessment recommendations
echo "DEBUG: Generating quality assessment recommendations..."

# Read key metrics for recommendations
MAPPING_RATE_NUM=$(grep "Mapping_rate_percent" "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt" 2>/dev/null | cut -f2 || echo "0")
MONO_PERCENTAGE=$(grep "Mono_nucleosome_percentage" "$FRAGMENT_LENGTH_STATS" | cut -f2)
MITO_PERC=$(grep "Mitochondrial_percentage" "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt" | cut -f2)
HIGH_QUALITY_BARCODES=$BARCODES_GT_500

cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_quality_recommendations.txt" << EOF
# Quality Assessment and Recommendations for $SAMPLE
# Generated: $(date)

=== OVERALL ASSESSMENT ===
Sample: $SAMPLE
Analysis Date: $(date)

=== KEY METRICS ===
Mapping Rate: ${MAPPING_RATE_NUM}%
Mono-nucleosome Percentage: ${MONO_PERCENTAGE}%
Mitochondrial Percentage: ${MITO_PERC}%
High-quality Barcodes (>500 frags): $HIGH_QUALITY_BARCODES

=== QUALITY RECOMMENDATIONS ===

MAPPING RATE:
$(if (( $(echo "$MAPPING_RATE_NUM > 60" | bc -l 2>/dev/null || echo 0) )); then
    echo "✅ GOOD: Mapping rate > 60% indicates good alignment quality"
elif (( $(echo "$MAPPING_RATE_NUM > 40" | bc -l 2>/dev/null || echo 0) )); then
    echo "⚠️  FAIR: Mapping rate 40-60% is acceptable but could be improved"
else
    echo "❌ POOR: Mapping rate < 40% suggests alignment issues"
fi)

FRAGMENT LENGTH DISTRIBUTION:
$(if (( $(echo "$MONO_PERCENTAGE > 40" | bc -l 2>/dev/null || echo 0) )); then
    echo "✅ GOOD: High mono-nucleosome percentage indicates good nucleosome structure"
else
    echo "⚠️  CHECK: Low mono-nucleosome percentage - verify fragment size selection"
fi)

MITOCHONDRIAL CONTENT:
$(if (( $(echo "$MITO_PERC < 20" | bc -l 2>/dev/null || echo 0) )); then
    echo "✅ GOOD: Low mitochondrial percentage indicates good nuclear enrichment"
elif (( $(echo "$MITO_PERC < 30" | bc -l 2>/dev/null || echo 0) )); then
    echo "⚠️  FAIR: Moderate mitochondrial content - still usable"
else
    echo "❌ POOR: High mitochondrial content suggests poor nuclear enrichment"
fi)

CELL QUALITY:
$(if [[ $HIGH_QUALITY_BARCODES -gt 1000 ]]; then
    echo "✅ GOOD: >1000 high-quality cells detected"
elif [[ $HIGH_QUALITY_BARCODES -gt 500 ]]; then
    echo "⚠️  FAIR: 500-1000 high-quality cells detected"
else
    echo "❌ POOR: <500 high-quality cells detected"
fi)

=== RECOMMENDATION FOR PEAK CALLING ===
$(
overall_score=0
if (( $(echo "$MAPPING_RATE_NUM > 50" | bc -l 2>/dev/null || echo 0) )); then ((overall_score++)); fi
if (( $(echo "$MONO_PERCENTAGE > 35" | bc -l 2>/dev/null || echo 0) )); then ((overall_score++)); fi
if (( $(echo "$MITO_PERC < 25" | bc -l 2>/dev/null || echo 0) )); then ((overall_score++)); fi
if [[ $HIGH_QUALITY_BARCODES -gt 500 ]]; then ((overall_score++)); fi

if [[ $overall_score -ge 3 ]]; then
    echo "✅ PROCEED: Data quality is sufficient for peak calling"
    echo "   Recommended to continue with step 5 (05_call_peaks.sh)"
elif [[ $overall_score -ge 2 ]]; then
    echo "⚠️  PROCEED WITH CAUTION: Data quality is marginal"
    echo "   Consider reviewing alignment parameters or filtering thresholds"
    echo "   You may still proceed with peak calling but monitor results carefully"
else
    echo "❌ REVIEW REQUIRED: Data quality issues detected"
    echo "   Recommend reviewing alignment parameters before peak calling"
    echo "   Consider re-processing with adjusted settings"
fi
)

=== FILES GENERATED ===
- Comprehensive QC summary: ${SAMPLE}_alignment_qc_summary.txt
- Fragment statistics: ${SAMPLE}_fragment_distribution.txt
- Length statistics: ${SAMPLE}_fragment_length_stats.txt
- Chromosome distribution: ${SAMPLE}_chromosome_distribution.txt
- Barcode quality: ${SAMPLE}_barcode_quality.txt
- QC plots: ${SAMPLE}_*.png (if R available)
- This recommendation file: ${SAMPLE}_quality_recommendations.txt
EOF

echo "Quality Assessment Summary:"
echo "=========================="
echo "Sample: $SAMPLE"
echo "Total Fragments: $(printf "%'d" $TOTAL_FRAGMENTS)"
echo "Unique Barcodes: $(printf "%'d" $UNIQUE_BARCODES)"
echo "Mapping Rate: ${MAPPING_RATE_NUM}%"
echo "Mitochondrial Content: ${MITO_PERC}%"
echo "High-Quality Cells (>500 frags): $(printf "%'d" $HIGH_QUALITY_BARCODES)"
echo ""
echo "📊 Detailed QC report: $OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_qc_summary.txt"
echo "📋 Quality recommendations: $OUTPUT_DIR/alignment_qc/${SAMPLE}_quality_recommendations.txt"

echo "Output files created in: $OUTPUT_DIR/alignment_qc/"
echo "  - Alignment QC summary: ${SAMPLE}_alignment_qc_summary.txt"
echo "  - Quality recommendations: ${SAMPLE}_quality_recommendations.txt"
echo "  - Fragment statistics: ${SAMPLE}_fragment_distribution.txt"
echo "  - Fragment length stats: ${SAMPLE}_fragment_length_stats.txt"
echo "  - Chromosome distribution: ${SAMPLE}_chromosome_distribution.txt"
echo "  - R plotting script: ${SAMPLE}_plot_alignment_qc.R"

echo "========================================="
echo "Step 4b complete for $SAMPLE"
echo "Alignment quality control analysis finished"
echo "End time: $(date)"
echo "========================================="