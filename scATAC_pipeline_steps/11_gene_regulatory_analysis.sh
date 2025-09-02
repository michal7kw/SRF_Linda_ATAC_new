#!/bin/bash
#SBATCH --job-name=grn_analysis
#SBATCH --output=logs/grn_analysis_%a.out
#SBATCH --error=logs/grn_analysis_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
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
echo "Step 11: Gene Regulatory Network Analysis for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create output directory
mkdir -p "$OUTPUT_DIR/grn_analysis"

# Check prerequisites
INTEGRATION_DIR="$OUTPUT_DIR/integration_prep"
PEAK_ANNOTATIONS="$INTEGRATION_DIR/${SAMPLE}_peak_annotations.tsv"
GENE_ACTIVITY="$INTEGRATION_DIR/${SAMPLE}_gene_activity.mtx"
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

if [[ ! -f "$PEAK_ANNOTATIONS" ]]; then
    echo "ERROR: Peak annotations not found: $PEAK_ANNOTATIONS"
    echo "Please run step 10 (10_integration_prep.sh) first"
    exit 1
fi

echo "DEBUG: Prerequisites verified"

# Create comprehensive R script for GRN analysis
cat > "$OUTPUT_DIR/grn_analysis/${SAMPLE}_grn_analysis.R" << 'EOF'
#!/usr/bin/env Rscript
# Gene Regulatory Network Analysis for scATAC-seq data
# Identifies TF binding sites, motif enrichment, and regulatory networks

suppressPackageStartupMessages({
    library(Matrix)
    library(GenomicRanges)
    library(motifmatchr)
    library(TFBSTools)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(chromVAR)
    library(SummarizedExperiment)
    library(ggplot2)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
    
    # Try to load JASPAR motifs
    if(!require(JASPAR2020, quietly = TRUE)) {
        cat("Installing JASPAR2020...\n")
        BiocManager::install("JASPAR2020", ask = FALSE)
        library(JASPAR2020)
    }
    
    if(!require(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)) {
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", ask = FALSE)
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    }
    
    if(!require(org.Mm.eg.db, quietly = TRUE)) {
        BiocManager::install("org.Mm.eg.db", ask = FALSE)
        library(org.Mm.eg.db)
    }
})

# Get sample name
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
    sample_name <- Sys.getenv("SAMPLE", "R26-Nestin-Ctrl-adult")
} else {
    sample_name <- args[1]
}

cat("Processing GRN analysis for:", sample_name, "\n")

# Set paths
output_dir <- "grn_analysis"
integration_dir <- "integration_prep"

# 1. Load peak data and annotations
cat("Loading peak data and annotations...\n")

# Load peak-barcode matrix
peak_matrix <- readMM(paste0(sample_name, "_peak_bc_matrix/matrix.mtx.gz"))
peak_features <- read.table(paste0(sample_name, "_peak_bc_matrix/features.tsv.gz"), 
                           header = FALSE, stringsAsFactors = FALSE)
peak_barcodes <- read.table(paste0(sample_name, "_peak_bc_matrix/barcodes.tsv.gz"), 
                           header = FALSE, stringsAsFactors = FALSE)

rownames(peak_matrix) <- peak_features[,1]
colnames(peak_matrix) <- peak_barcodes[,1]

# Load peak annotations
peak_annotations <- read.table(file.path(integration_dir, paste0(sample_name, "_peak_annotations.tsv")),
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Loaded", nrow(peak_matrix), "peaks x", ncol(peak_matrix), "cells\n")

# 2. Create GenomicRanges object for peaks
cat("Creating genomic ranges for peaks...\n")
peak_coords <- rownames(peak_matrix)
peak_chr <- gsub(":.*", "", peak_coords)
peak_ranges <- gsub(".*:", "", peak_coords)
peak_start <- as.numeric(gsub("-.*", "", peak_ranges))
peak_end <- as.numeric(gsub(".*-", "", peak_ranges))

peak_gr <- GRanges(
    seqnames = peak_chr,
    ranges = IRanges(start = peak_start, end = peak_end),
    peak_id = peak_coords
)

# Filter for standard chromosomes only
standard_chroms <- paste0("chr", c(1:19, "X", "Y"))
peak_gr <- peak_gr[seqnames(peak_gr) %in% standard_chroms]
peak_matrix <- peak_matrix[peak_gr$peak_id, ]

cat("Filtered to", length(peak_gr), "peaks on standard chromosomes\n")

# 3. Get motif database
cat("Loading JASPAR motif database...\n")
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"
opts[["matrixtype"]] <- "PFM"
motifs <- getMatrixSet(JASPAR2020, opts)

# Filter for mouse TFs
mouse_tfs <- c("Ascl1", "Dlx2", "Emx2", "Foxg1", "Gsx2", "Hes1", "Hes5", 
               "Hey1", "Hey2", "Lhx2", "Mash1", "Neurog2", "Olig2", "Pax6", 
               "Sox2", "Tbr1", "Tbr2", "Tlx", "Zic1", "Zic3")

# Keep motifs for relevant TFs (expand this list based on your specific interests)
relevant_motifs <- motifs[TFBSTools::name(motifs) %in% mouse_tfs]
cat("Using", length(relevant_motifs), "motifs for analysis\n")

# 4. Find motif matches in peaks
cat("Finding motif matches in accessible peaks...\n")
genome <- BSgenome.Mmusculus.UCSC.mm10

# This step can be memory intensive, so we'll process in chunks if needed
chunk_size <- min(5000, length(peak_gr))
n_chunks <- ceiling(length(peak_gr) / chunk_size)

motif_matches_list <- list()
for(i in 1:n_chunks) {
    cat("Processing chunk", i, "of", n_chunks, "\n")
    
    start_idx <- ((i-1) * chunk_size) + 1
    end_idx <- min(i * chunk_size, length(peak_gr))
    
    peak_chunk <- peak_gr[start_idx:end_idx]
    
    # Get sequences
    peak_seqs <- getSeq(genome, peak_chunk)
    
    # Find motif matches
    motif_matches_chunk <- matchMotifs(relevant_motifs, peak_seqs, 
                                      genome = genome, p.cutoff = 1e-4)
    
    motif_matches_list[[i]] <- assay(motif_matches_chunk)
}

# Combine results
motif_matches <- do.call(rbind, motif_matches_list)
colnames(motif_matches) <- TFBSTools::name(relevant_motifs)
rownames(motif_matches) <- peak_gr$peak_id

cat("Found motif matches in", nrow(motif_matches), "peaks for", ncol(motif_matches), "TFs\n")

# 5. Calculate TF activity scores using chromVAR approach
cat("Calculating TF activity scores...\n")

# Create SummarizedExperiment object
peak_se <- SummarizedExperiment(
    assays = list(counts = peak_matrix[rownames(motif_matches), ]),
    rowRanges = peak_gr[match(rownames(motif_matches), peak_gr$peak_id)],
    colData = DataFrame(cell_id = colnames(peak_matrix))
)

# Add motif matches
rowData(peak_se)$motif_matches <- motif_matches

# Calculate background peaks
bg_peaks <- getBackgroundPeaks(peak_se, niterations = 200)

# Calculate deviations (TF activity)
dev_scores <- computeDeviations(object = peak_se, 
                               annotations = motif_matches,
                               background_peaks = bg_peaks)

# Extract z-scores (TF activity per cell)
tf_activity <- assays(dev_scores)[["z"]]
rownames(tf_activity) <- colnames(motif_matches)

cat("Calculated TF activity scores:", nrow(tf_activity), "TFs x", ncol(tf_activity), "cells\n")

# 6. Identify differentially active TFs (if comparing conditions)
cat("Analyzing TF activity patterns...\n")

# Calculate TF activity statistics
tf_stats <- data.frame(
    TF = rownames(tf_activity),
    mean_activity = rowMeans(tf_activity),
    sd_activity = apply(tf_activity, 1, sd),
    max_activity = apply(tf_activity, 1, max),
    min_activity = apply(tf_activity, 1, min)
)

tf_stats$cv <- tf_stats$sd_activity / abs(tf_stats$mean_activity)
tf_stats <- tf_stats[order(tf_stats$cv, decreasing = TRUE), ]

# 7. Peak-to-gene linking for regulatory networks
cat("Creating peak-to-gene regulatory links...\n")

# Use existing peak annotations
peak_gene_links <- peak_annotations[peak_annotations$associated_genes != "none", ]

# Create regulatory network edges
regulatory_edges <- data.frame()

for(i in 1:nrow(peak_gene_links)) {
    peak_id <- peak_gene_links$peak_id[i]
    genes <- unlist(strsplit(peak_gene_links$associated_genes[i], ","))
    
    # Find TFs with motifs in this peak
    peak_motifs <- names(which(motif_matches[peak_id, ]))
    
    if(length(peak_motifs) > 0 && length(genes) > 0) {
        # Create TF -> target gene edges
        for(tf in peak_motifs) {
            for(gene in genes) {
                regulatory_edges <- rbind(regulatory_edges, data.frame(
                    TF = tf,
                    target_gene = gene,
                    peak_id = peak_id,
                    regulation_type = peak_gene_links$association_types[i],
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
}

cat("Identified", nrow(regulatory_edges), "potential regulatory interactions\n")

# 8. Save results
cat("Saving results...\n")

# Save TF activity matrix
write.table(tf_activity, 
           file.path(output_dir, paste0(sample_name, "_tf_activity.tsv")),
           sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Save TF statistics
write.table(tf_stats,
           file.path(output_dir, paste0(sample_name, "_tf_statistics.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save regulatory network
write.table(regulatory_edges,
           file.path(output_dir, paste0(sample_name, "_regulatory_network.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save motif matches
motif_matches_df <- data.frame(
    peak_id = rep(rownames(motif_matches), ncol(motif_matches)),
    TF = rep(colnames(motif_matches), each = nrow(motif_matches)),
    has_motif = as.vector(motif_matches)
)
motif_matches_df <- motif_matches_df[motif_matches_df$has_motif, ]

write.table(motif_matches_df[, c("peak_id", "TF")],
           file.path(output_dir, paste0(sample_name, "_motif_matches.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 9. Generate visualizations
cat("Generating visualizations...\n")

# TF activity heatmap (top variable TFs)
top_tfs <- head(tf_stats$TF, 20)
tf_activity_subset <- tf_activity[top_tfs, ]

# Sample cells for visualization if too many
if(ncol(tf_activity_subset) > 1000) {
    sampled_cells <- sample(ncol(tf_activity_subset), 1000)
    tf_activity_subset <- tf_activity_subset[, sampled_cells]
}

png(file.path(output_dir, paste0(sample_name, "_tf_activity_heatmap.png")), 
    width = 12, height = 10, units = "in", res = 300)

# Create heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(tf_activity_subset, 
        name = "TF Activity\n(z-score)",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 12)))

dev.off()

# TF activity distribution
tf_stats_plot <- tf_stats %>% head(15)
p1 <- ggplot(tf_stats_plot, aes(x = reorder(TF, cv), y = cv)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    coord_flip() +
    labs(title = paste("Most Variable TFs -", sample_name),
         x = "Transcription Factor", 
         y = "Coefficient of Variation") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_tf_variability.png")), 
       p1, width = 10, height = 8, dpi = 300)

# Regulatory network summary
network_summary <- regulatory_edges %>%
    group_by(TF) %>%
    summarise(n_targets = n_distinct(target_gene)) %>%
    arrange(desc(n_targets)) %>%
    head(15)

p2 <- ggplot(network_summary, aes(x = reorder(TF, n_targets), y = n_targets)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    coord_flip() +
    labs(title = paste("TFs with Most Target Genes -", sample_name),
         x = "Transcription Factor", 
         y = "Number of Target Genes") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_tf_targets.png")), 
       p2, width = 10, height = 8, dpi = 300)

# 10. Create summary report
cat("Creating summary report...\n")

grn_summary <- data.frame(
    metric = c("total_peaks_analyzed", "total_cells", "motifs_searched", 
               "peaks_with_motifs", "tfs_with_activity", "regulatory_edges",
               "mean_tf_activity", "most_active_tf", "most_variable_tf"),
    value = c(nrow(motif_matches), ncol(tf_activity), ncol(motif_matches),
              sum(rowSums(motif_matches) > 0), nrow(tf_activity), nrow(regulatory_edges),
              round(mean(abs(tf_stats$mean_activity)), 3),
              tf_stats$TF[which.max(abs(tf_stats$mean_activity))],
              tf_stats$TF[1])
)

write.table(grn_summary,
           file.path(output_dir, paste0(sample_name, "_grn_summary.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("GRN analysis complete!\n")
cat("Files generated:\n")
cat("  - TF activity matrix:", paste0(sample_name, "_tf_activity.tsv"), "\n")
cat("  - TF statistics:", paste0(sample_name, "_tf_statistics.tsv"), "\n")
cat("  - Regulatory network:", paste0(sample_name, "_regulatory_network.tsv"), "\n")
cat("  - Motif matches:", paste0(sample_name, "_motif_matches.tsv"), "\n")
cat("  - Summary report:", paste0(sample_name, "_grn_summary.tsv"), "\n")
cat("  - Visualizations: heatmaps and bar plots\n")

# Print key findings
cat("\nKey Findings:\n")
cat("=============\n")
cat("Most variable TFs:\n")
print(head(tf_stats[, c("TF", "cv")], 5))
cat("\nTop regulatory hubs:\n")
print(head(network_summary, 5))
EOF

# Make R script executable
chmod +x "$OUTPUT_DIR/grn_analysis/${SAMPLE}_grn_analysis.R"

# Create a comparative analysis script for control vs mutant
if [[ "$SAMPLE" == "R26-Nestin-Ctrl-adult" ]]; then
cat > "$OUTPUT_DIR/grn_analysis/comparative_grn_analysis.R" << 'EOF'
#!/usr/bin/env Rscript
# Comparative GRN Analysis: Control vs Mutant
# Identifies differentially active TFs and disrupted regulatory networks

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
})

cat("Running comparative GRN analysis...\n")

ctrl_sample <- "R26-Nestin-Ctrl-adult"
mut_sample <- "R26-Nestin-Mut-adult"
output_dir <- "grn_analysis"

# Load TF activity data
ctrl_tf <- tryCatch({
    read.table(file.path(output_dir, paste0(ctrl_sample, "_tf_activity.tsv")),
               header = TRUE, sep = "\t", row.names = 1)
}, error = function(e) {
    cat("Control data not found\n")
    return(NULL)
})

mut_tf <- tryCatch({
    read.table(file.path(output_dir, paste0(mut_sample, "_tf_activity.tsv")),
               header = TRUE, sep = "\t", row.names = 1)
}, error = function(e) {
    cat("Mutant data not found\n")
    return(NULL)
})

if(is.null(ctrl_tf) || is.null(mut_tf)) {
    stop("Cannot perform comparative analysis - missing data files")
}

# Find common TFs
common_tfs <- intersect(rownames(ctrl_tf), rownames(mut_tf))
cat("Comparing", length(common_tfs), "common TFs\n")

if(length(common_tfs) == 0) {
    stop("No common TFs found between samples")
}

# Calculate mean TF activity per condition
ctrl_means <- rowMeans(ctrl_tf[common_tfs, ])
mut_means <- rowMeans(mut_tf[common_tfs, ])

# Differential TF activity
tf_comparison <- data.frame(
    TF = common_tfs,
    ctrl_activity = ctrl_means,
    mut_activity = mut_means,
    log2fc = log2((mut_means + 1) / (ctrl_means + 1)),
    difference = mut_means - ctrl_means
)

tf_comparison$abs_log2fc <- abs(tf_comparison$log2fc)
tf_comparison <- tf_comparison[order(tf_comparison$abs_log2fc, decreasing = TRUE), ]

# Identify significantly changed TFs
tf_comparison$direction <- ifelse(tf_comparison$log2fc > 0.5, "Up in Mutant",
                                ifelse(tf_comparison$log2fc < -0.5, "Down in Mutant", "Unchanged"))

cat("TF activity changes:\n")
cat("Up in Mutant:", sum(tf_comparison$direction == "Up in Mutant"), "\n")
cat("Down in Mutant:", sum(tf_comparison$direction == "Down in Mutant"), "\n")
cat("Unchanged:", sum(tf_comparison$direction == "Unchanged"), "\n")

# Save comparative results
write.table(tf_comparison,
           file.path(output_dir, "ctrl_vs_mut_tf_comparison.tsv"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Generate comparative plots
# Volcano plot
p1 <- ggplot(tf_comparison, aes(x = log2fc, y = abs_log2fc, color = direction)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(aes(label = ifelse(abs_log2fc > 1, TF, "")), 
              vjust = -0.5, size = 3) +
    scale_color_manual(values = c("Up in Mutant" = "red", 
                                 "Down in Mutant" = "blue", 
                                 "Unchanged" = "gray")) +
    labs(title = "Differential TF Activity: Control vs Mutant",
         x = "Log2 Fold Change (Mutant/Control)",
         y = "Absolute Log2 Fold Change",
         color = "Direction") +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(file.path(output_dir, "tf_activity_volcano_plot.png"), 
       p1, width = 10, height = 8, dpi = 300)

# Correlation plot
p2 <- ggplot(tf_comparison, aes(x = ctrl_activity, y = mut_activity)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_text(aes(label = ifelse(abs_log2fc > 1, TF, "")), 
              vjust = -0.5, hjust = 0.5, size = 3) +
    labs(title = "TF Activity Correlation: Control vs Mutant",
         x = "Control TF Activity (mean z-score)",
         y = "Mutant TF Activity (mean z-score)") +
    theme_minimal()

ggsave(file.path(output_dir, "tf_activity_correlation.png"), 
       p2, width = 8, height = 8, dpi = 300)

cat("Comparative analysis complete!\n")
cat("Top differentially active TFs:\n")
print(head(tf_comparison[, c("TF", "log2fc", "direction")], 10))
EOF

chmod +x "$OUTPUT_DIR/grn_analysis/comparative_grn_analysis.R"
fi

# Run the GRN analysis
echo "DEBUG: Running gene regulatory network analysis..."
cd "$OUTPUT_DIR"
export SAMPLE="$SAMPLE"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R"
    exit 1
fi

# Run the main GRN analysis script
echo "DEBUG: Starting GRN analysis for $SAMPLE..."
Rscript "grn_analysis/${SAMPLE}_grn_analysis.R" "$SAMPLE"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: GRN analysis completed successfully"
    
    # Display summary if available
    SUMMARY_FILE="grn_analysis/${SAMPLE}_grn_summary.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "GRN Analysis Summary:"
        cat "$SUMMARY_FILE"
    fi
else
    echo "WARNING: GRN analysis failed. Check R dependencies."
    echo "Required packages: motifmatchr, TFBSTools, BSgenome.Mmusculus.UCSC.mm10, chromVAR"
    echo "Install with BiocManager::install(c('motifmatchr', 'TFBSTools', 'BSgenome.Mmusculus.UCSC.mm10', 'chromVAR'))"
fi

echo "Output files created in: $OUTPUT_DIR/grn_analysis/"
echo "  - TF activity scores: ${SAMPLE}_tf_activity.tsv"
echo "  - TF statistics: ${SAMPLE}_tf_statistics.tsv"
echo "  - Regulatory network: ${SAMPLE}_regulatory_network.tsv"
echo "  - Motif matches: ${SAMPLE}_motif_matches.tsv"
echo "  - Analysis summary: ${SAMPLE}_grn_summary.tsv"
echo "  - Visualizations: heatmaps and plots"

echo "========================================="
echo "Step 11 complete for $SAMPLE"
echo "Gene regulatory network analysis finished"
echo "End time: $(date)"
echo "========================================="