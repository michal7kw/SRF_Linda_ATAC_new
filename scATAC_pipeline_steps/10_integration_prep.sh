#!/bin/bash
#SBATCH --job-name=integration_prep
#SBATCH --output=logs/integration_prep_%a.out
#SBATCH --error=logs/integration_prep_%a.err
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
echo "Step 10: Integration Preparation for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create output directory
mkdir -p "$OUTPUT_DIR/integration_prep"

# Check prerequisites
PEAK_BC_MATRIX="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

if [[ ! -d "$PEAK_BC_MATRIX" ]]; then
    echo "ERROR: Peak-barcode matrix directory not found: $PEAK_BC_MATRIX"
    echo "Please run step 6 (06_peak_cell_matrix.sh) first"
    exit 1
fi

if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    exit 1
fi

echo "DEBUG: Prerequisites verified"

# Create comprehensive R script for integration preparation
cat > "$OUTPUT_DIR/integration_prep/${SAMPLE}_integration_prep.R" << 'EOF'
#!/usr/bin/env Rscript
# Integration Preparation Script for scATAC-seq data
# Prepares data for integration with scRNA-seq and GRN analysis

suppressPackageStartupMessages({
    library(Matrix)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(dplyr)
    library(ggplot2)
    if(!require(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)) {
        cat("Installing TxDb.Mmusculus.UCSC.mm10.knownGene...\n")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", ask = FALSE)
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    }
    if(!require(org.Mm.eg.db, quietly = TRUE)) {
        cat("Installing org.Mm.eg.db...\n")
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

cat("Processing sample:", sample_name, "\n")

# Set paths
matrix_dir <- paste0(sample_name, "_peak_bc_matrix")
output_dir <- "integration_prep"

# 1. Load peak-barcode matrix
cat("Loading peak-barcode matrix...\n")
matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")
features_path <- file.path(matrix_dir, "features.tsv.gz")
barcodes_path <- file.path(matrix_dir, "barcodes.tsv.gz")

peak_matrix <- readMM(matrix_path)
features <- read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(barcodes_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

rownames(peak_matrix) <- features[,1]
colnames(peak_matrix) <- barcodes[,1]

cat("Matrix loaded:", nrow(peak_matrix), "peaks x", ncol(peak_matrix), "cells\n")

# 2. Create gene activity scores
cat("Creating gene activity scores for integration...\n")

# Load gene annotations
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)

# Parse peak coordinates
peak_coords <- rownames(peak_matrix)
peak_chr <- gsub(":.*", "", peak_coords)
peak_ranges <- gsub(".*:", "", peak_coords)
peak_start <- as.numeric(gsub("-.*", "", peak_ranges))
peak_end <- as.numeric(gsub(".*-", "", peak_ranges))

# Create GenomicRanges object for peaks
peak_gr <- GRanges(
    seqnames = peak_chr,
    ranges = IRanges(start = peak_start, end = peak_end)
)
names(peak_gr) <- peak_coords

cat("Created", length(peak_gr), "peak ranges\n")

# Find overlaps between peaks and gene bodies + promoters
cat("Finding peak-gene associations...\n")

# Gene bodies
gene_overlaps <- findOverlaps(peak_gr, genes)
gene_body_associations <- data.frame(
    peak = names(peak_gr)[queryHits(gene_overlaps)],
    gene_id = names(genes)[subjectHits(gene_overlaps)],
    type = "gene_body"
)

# Promoter regions (TSS +/- 2kb)
promoters <- promoters(genes, upstream = 2000, downstream = 2000)
promoter_overlaps <- findOverlaps(peak_gr, promoters)
promoter_associations <- data.frame(
    peak = names(peak_gr)[queryHits(promoter_overlaps)],
    gene_id = names(promoters)[subjectHits(promoter_overlaps)],
    type = "promoter"
)

# Combine associations
all_associations <- rbind(gene_body_associations, promoter_associations)
cat("Found", nrow(all_associations), "peak-gene associations\n")

# Create gene activity matrix
cat("Computing gene activity scores...\n")
unique_genes <- unique(all_associations$gene_id)
gene_activity_matrix <- Matrix(0, nrow = length(unique_genes), 
                              ncol = ncol(peak_matrix), sparse = TRUE)
rownames(gene_activity_matrix) <- unique_genes
colnames(gene_activity_matrix) <- colnames(peak_matrix)

# Sum peak accessibility for each gene
for(gene in unique_genes) {
    associated_peaks <- all_associations$peak[all_associations$gene_id == gene]
    if(length(associated_peaks) > 0) {
        # Find peaks that exist in our matrix
        valid_peaks <- intersect(associated_peaks, rownames(peak_matrix))
        if(length(valid_peaks) > 0) {
            if(length(valid_peaks) == 1) {
                gene_activity_matrix[gene, ] <- peak_matrix[valid_peaks, ]
            } else {
                gene_activity_matrix[gene, ] <- Matrix::colSums(peak_matrix[valid_peaks, ])
            }
        }
    }
}

# Remove genes with no activity
gene_sums <- Matrix::rowSums(gene_activity_matrix)
active_genes <- gene_sums > 0
gene_activity_matrix <- gene_activity_matrix[active_genes, ]

cat("Gene activity matrix:", nrow(gene_activity_matrix), "genes x", ncol(gene_activity_matrix), "cells\n")

# 3. Convert gene IDs to symbols
cat("Converting gene IDs to symbols...\n")
gene_symbols <- tryCatch({
    mapIds(org.Mm.eg.db, keys = rownames(gene_activity_matrix), 
           column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
}, error = function(e) {
    cat("Warning: Could not convert gene IDs to symbols\n")
    setNames(rownames(gene_activity_matrix), rownames(gene_activity_matrix))
})

# Remove NAs and duplicates
valid_symbols <- !is.na(gene_symbols) & !duplicated(gene_symbols)
gene_activity_matrix <- gene_activity_matrix[valid_symbols, ]
rownames(gene_activity_matrix) <- gene_symbols[valid_symbols]

cat("Final gene activity matrix:", nrow(gene_activity_matrix), "genes x", ncol(gene_activity_matrix), "cells\n")

# 4. Save gene activity matrix in multiple formats
cat("Saving gene activity matrix...\n")

# Save as Matrix Market format
Matrix::writeMM(gene_activity_matrix, 
               file.path(output_dir, paste0(sample_name, "_gene_activity.mtx")))

# Save gene names
write.table(data.frame(gene_symbol = rownames(gene_activity_matrix)),
           file.path(output_dir, paste0(sample_name, "_gene_activity_features.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save cell barcodes (same as peak matrix)
write.table(data.frame(barcode = colnames(gene_activity_matrix)),
           file.path(output_dir, paste0(sample_name, "_gene_activity_barcodes.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 5. Create Seurat-compatible objects
cat("Creating Seurat-compatible metadata...\n")

# Calculate per-cell metrics
cell_metrics <- data.frame(
    barcode = colnames(peak_matrix),
    n_peaks = Matrix::colSums(peak_matrix > 0),
    total_fragments = Matrix::colSums(peak_matrix),
    n_genes = Matrix::colSums(gene_activity_matrix > 0),
    total_gene_activity = Matrix::colSums(gene_activity_matrix)
)

# Calculate additional QC metrics
cell_metrics$log10_total_fragments <- log10(cell_metrics$total_fragments + 1)
cell_metrics$peaks_per_fragment <- cell_metrics$n_peaks / cell_metrics$total_fragments

# Save cell metadata
write.table(cell_metrics,
           file.path(output_dir, paste0(sample_name, "_cell_metadata.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 6. Create peak annotations for regulatory analysis
cat("Creating peak annotations...\n")
peak_annotations <- data.frame(
    peak_id = names(peak_gr),
    chr = as.character(seqnames(peak_gr)),
    start = start(peak_gr),
    end = end(peak_gr),
    width = width(peak_gr)
)

# Add gene associations
peak_gene_map <- all_associations %>%
    group_by(peak) %>%
    summarise(
        associated_genes = paste(gene_id, collapse = ","),
        n_genes = n(),
        association_types = paste(unique(type), collapse = ",")
    )

peak_annotations <- merge(peak_annotations, peak_gene_map, 
                         by.x = "peak_id", by.y = "peak", all.x = TRUE)

# Fill NAs
peak_annotations$associated_genes[is.na(peak_annotations$associated_genes)] <- "none"
peak_annotations$n_genes[is.na(peak_annotations$n_genes)] <- 0
peak_annotations$association_types[is.na(peak_annotations$association_types)] <- "intergenic"

write.table(peak_annotations,
           file.path(output_dir, paste0(sample_name, "_peak_annotations.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 7. Create integration summary
cat("Creating integration summary...\n")
integration_summary <- data.frame(
    metric = c("total_peaks", "total_cells", "active_genes", "mean_peaks_per_cell",
               "mean_fragments_per_cell", "mean_gene_activity_per_cell",
               "peaks_with_gene_associations", "genes_with_peak_associations"),
    value = c(nrow(peak_matrix), ncol(peak_matrix), nrow(gene_activity_matrix),
              round(mean(cell_metrics$n_peaks), 2),
              round(mean(cell_metrics$total_fragments), 2),
              round(mean(cell_metrics$total_gene_activity), 2),
              sum(peak_annotations$n_genes > 0),
              length(unique(unlist(strsplit(peak_annotations$associated_genes[
                peak_annotations$associated_genes != "none"], ",")))))
)

write.table(integration_summary,
           file.path(output_dir, paste0(sample_name, "_integration_summary.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 8. Generate QC plots for integration
cat("Generating integration QC plots...\n")

# Peaks vs gene activity per cell
p1 <- ggplot(cell_metrics, aes(x = n_peaks, y = n_genes)) +
    geom_point(alpha = 0.6, size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    labs(title = paste("Peaks vs Gene Activity -", sample_name),
         x = "Number of Accessible Peaks", 
         y = "Number of Active Genes") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_peaks_vs_genes.png")), 
       p1, width = 8, height = 6, dpi = 300)

# Distribution of gene activity scores
gene_activity_dist <- data.frame(
    gene = rownames(gene_activity_matrix),
    total_activity = Matrix::rowSums(gene_activity_matrix),
    n_cells_active = Matrix::rowSums(gene_activity_matrix > 0)
)

p2 <- ggplot(gene_activity_dist, aes(x = n_cells_active, y = total_activity)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste("Gene Activity Distribution -", sample_name),
         x = "Number of Cells with Activity (log10)",
         y = "Total Activity Score (log10)") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_gene_activity_dist.png")), 
       p2, width = 8, height = 6, dpi = 300)

cat("Integration preparation complete!\n")
cat("Files generated:\n")
cat("  - Gene activity matrix: ", paste0(sample_name, "_gene_activity.mtx"), "\n")
cat("  - Gene features: ", paste0(sample_name, "_gene_activity_features.tsv"), "\n")
cat("  - Cell barcodes: ", paste0(sample_name, "_gene_activity_barcodes.tsv"), "\n")
cat("  - Cell metadata: ", paste0(sample_name, "_cell_metadata.tsv"), "\n")
cat("  - Peak annotations: ", paste0(sample_name, "_peak_annotations.tsv"), "\n")
cat("  - Integration summary: ", paste0(sample_name, "_integration_summary.tsv"), "\n")
cat("  - QC plots: peaks vs genes, gene activity distribution\n")
EOF

# Make R script executable
chmod +x "$OUTPUT_DIR/integration_prep/${SAMPLE}_integration_prep.R"

# Create a simple Seurat integration helper script
cat > "$OUTPUT_DIR/integration_prep/${SAMPLE}_load_seurat.R" << 'EOF'
#!/usr/bin/env Rscript
# Helper script to load scATAC data into Seurat for integration

# This script creates Seurat objects from the processed scATAC data
# Usage: Rscript load_seurat.R sample_name

suppressPackageStartupMessages({
    library(Seurat)
    library(Signac)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
    stop("Please provide sample name as argument")
}

sample_name <- args[1]
output_dir <- "integration_prep"

cat("Loading", sample_name, "data into Seurat...\n")

# Load peak matrix
peak_matrix <- readMM(paste0(sample_name, "_peak_bc_matrix/matrix.mtx.gz"))
peak_features <- read.table(paste0(sample_name, "_peak_bc_matrix/features.tsv.gz"), 
                           header = FALSE, stringsAsFactors = FALSE)
peak_barcodes <- read.table(paste0(sample_name, "_peak_bc_matrix/barcodes.tsv.gz"), 
                           header = FALSE, stringsAsFactors = FALSE)

rownames(peak_matrix) <- peak_features[,1]
colnames(peak_matrix) <- peak_barcodes[,1]

# Load gene activity matrix
gene_matrix <- readMM(file.path(output_dir, paste0(sample_name, "_gene_activity.mtx")))
gene_features <- read.table(file.path(output_dir, paste0(sample_name, "_gene_activity_features.tsv")),
                           header = FALSE, stringsAsFactors = FALSE)
gene_barcodes <- read.table(file.path(output_dir, paste0(sample_name, "_gene_activity_barcodes.tsv")),
                           header = FALSE, stringsAsFactors = FALSE)

rownames(gene_matrix) <- gene_features[,1]
colnames(gene_matrix) <- gene_barcodes[,1]

# Load metadata
metadata <- read.table(file.path(output_dir, paste0(sample_name, "_cell_metadata.tsv")),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$barcode

cat("Creating Seurat object...\n")

# Create ChromatinAssay for peaks
chrom_assay <- CreateChromatinAssay(
    counts = peak_matrix,
    sep = c(":", "-")
)

# Create Seurat object
atac_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
)

# Add gene activity as a separate assay
atac_obj[["RNA"]] <- CreateAssayObject(counts = gene_matrix)

# Set default assay
DefaultAssay(atac_obj) <- "peaks"

cat("Seurat object created successfully!\n")
cat("Assays:", names(atac_obj@assays), "\n")
cat("Cells:", ncol(atac_obj), "\n")
cat("Peaks:", nrow(atac_obj[["peaks"]]), "\n")
cat("Genes:", nrow(atac_obj[["RNA"]]), "\n")

# Save Seurat object
saveRDS(atac_obj, file.path(output_dir, paste0(sample_name, "_seurat_object.rds")))
cat("Seurat object saved to:", file.path(output_dir, paste0(sample_name, "_seurat_object.rds")), "\n")
EOF

chmod +x "$OUTPUT_DIR/integration_prep/${SAMPLE}_load_seurat.R"

# Run the integration preparation
echo "DEBUG: Running integration preparation analysis..."
cd "$OUTPUT_DIR"
export SAMPLE="$SAMPLE"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R"
    exit 1
fi

# Run the main integration prep script
echo "DEBUG: Starting integration preparation for $SAMPLE..."
Rscript "integration_prep/${SAMPLE}_integration_prep.R" "$SAMPLE"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: Integration preparation completed successfully"
    
    # Display summary if available
    SUMMARY_FILE="integration_prep/${SAMPLE}_integration_summary.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "Integration Summary:"
        cat "$SUMMARY_FILE"
    fi
else
    echo "WARNING: Integration preparation failed. Check R dependencies."
    echo "Required packages: Matrix, GenomicRanges, GenomicFeatures, rtracklayer"
    echo "Install with BiocManager::install(c('GenomicRanges', 'GenomicFeatures', 'rtracklayer'))"
fi

echo "Output files created in: $OUTPUT_DIR/integration_prep/"
echo "  - Gene activity matrix: ${SAMPLE}_gene_activity.mtx"
echo "  - Cell metadata: ${SAMPLE}_cell_metadata.tsv"
echo "  - Peak annotations: ${SAMPLE}_peak_annotations.tsv"
echo "  - Integration summary: ${SAMPLE}_integration_summary.tsv"
echo "  - Seurat helper script: ${SAMPLE}_load_seurat.R"
echo "  - QC plots: peaks vs genes, gene activity distribution"

echo "========================================="
echo "Step 10 complete for $SAMPLE"
echo "Integration preparation finished"
echo "End time: $(date)"
echo "========================================="