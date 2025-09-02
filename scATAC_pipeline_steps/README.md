# scATAC-seq Processing Pipeline

This directory contains a comprehensive single-cell ATAC-seq (scATAC-seq) processing pipeline for R26-Nestin Control vs Mutant analysis. The pipeline processes chromatin accessibility data to identify disrupted gene regulatory networks.

## Pipeline Overview

The pipeline consists of 12 scripts that should be executed in the following order:

```
01_extract_barcodes.sh ✅ (COMPLETED)
        ↓
02_test_barcodes.sh ✅ (COMPLETED)  
        ↓
03_validate_counts.sh ❌ (PENDING)
        ↓
04_chromap_alignment.sh ✅ (COMPLETED)
        ↓
05_call_peaks.sh ❌ (PENDING)
        ↓
06_peak_cell_matrix.sh ❌ (PENDING)
        ↓
┌─────────────────────────────────────┐
│     Advanced Analysis Scripts      │
│     (Can run in parallel)          │
├─────────────────────────────────────┤
│ 07_create_bigwig.sh                 │
│ 08_qc_metrics.sh                    │ 
│ 09_dimensionality_reduction.sh      │
│ 10_integration_prep.sh              │
│ 11_gene_regulatory_analysis.sh      │
└─────────────────────────────────────┘
        ↓
12_multiqc_report.sh (Generate final QC report)
```

## Script Descriptions

### Core Processing Scripts (Sequential)

#### 1. **01_extract_barcodes.sh** ✅ COMPLETED
- **Purpose**: Extract and validate cell barcodes from FASTQ files
- **Input**: Raw FASTQ files
- **Output**: Barcode lists in `qc/` directory
- **Dependencies**: None

#### 2. **02_test_barcodes.sh** ✅ COMPLETED  
- **Purpose**: Test barcode quality and create filtered barcode lists
- **Input**: Extracted barcodes from step 1
- **Output**: Cleaned barcode files in `qc/`
- **Dependencies**: 01_extract_barcodes.sh

#### 3. **03_validate_counts.sh** ❌ PENDING
- **Purpose**: Validate fragment counts and barcode statistics
- **Input**: Barcode files and fragment data
- **Output**: Count validation reports
- **Dependencies**: 02_test_barcodes.sh

#### 4. **04_chromap_alignment.sh** ✅ COMPLETED
- **Purpose**: Align reads using chromap and generate fragment files
- **Input**: FASTQ files and validated barcodes
- **Output**: 
  - Fragment files in `fragments/`
  - BED files with aligned reads
  - Chromap alignment logs
- **Dependencies**: 02_test_barcodes.sh (can run parallel to 03)

#### 5. **05_call_peaks.sh** ❌ PENDING
- **Purpose**: Call accessible chromatin peaks using MACS2
- **Input**: Fragment files or BED files from step 4
- **Output**: Peak files in `peaks/` directory
- **Dependencies**: 04_chromap_alignment.sh

#### 6. **06_peak_cell_matrix.sh** ❌ PENDING
- **Purpose**: Create peak-by-cell count matrix in 10X format
- **Input**: Peak files and fragment files
- **Output**: 
  - Peak-barcode matrix in `{sample}_peak_bc_matrix/`
  - `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`
- **Dependencies**: 05_call_peaks.sh

### Advanced Analysis Scripts (Can run in parallel after step 6)

#### 7. **07_create_bigwig.sh** ✨ NEW
- **Purpose**: Generate BigWig coverage tracks for visualization
- **Input**: BED files from step 4
- **Output**: 
  - Raw coverage BigWig files in `bigwig/`
  - RPM-normalized BigWig files
- **Dependencies**: 04_chromap_alignment.sh
- **Environment**: `conda activate bigwig`

#### 8. **08_qc_metrics.sh** ✨ NEW  
- **Purpose**: Comprehensive quality control analysis
- **Input**: Fragment files, peaks, and matrices
- **Output**: 
  - QC metrics in `qc_metrics/`
  - Fragment length distributions
  - TSS enrichment scores
  - Library complexity analysis
  - R plotting scripts
- **Dependencies**: 05_call_peaks.sh, 06_peak_cell_matrix.sh

#### 9. **09_dimensionality_reduction.sh** ✨ NEW
- **Purpose**: Perform LSI, PCA, UMAP, and clustering
- **Input**: Peak-barcode matrix from step 6
- **Output**: 
  - Dimensionality reduction results in `dimensionality_reduction/`
  - LSI, PCA, UMAP embeddings
  - Cluster assignments
  - Visualization plots
- **Dependencies**: 06_peak_cell_matrix.sh
- **Resources**: 16 CPUs, 64GB RAM, 8h

#### 10. **10_integration_prep.sh** ✨ NEW
- **Purpose**: Prepare data for scRNA-seq integration
- **Input**: Peak-barcode matrix and peak annotations
- **Output**: 
  - Gene activity scores in `integration_prep/`
  - Cell metadata for integration
  - Peak-to-gene annotations
  - Seurat-compatible objects
- **Dependencies**: 06_peak_cell_matrix.sh
- **Key Feature**: Creates gene activity scores for multi-modal integration

#### 11. **11_gene_regulatory_analysis.sh** ✨ NEW
- **Purpose**: Gene regulatory network (GRN) analysis and TF motif discovery
- **Input**: Peak matrices and integration data
- **Output**: 
  - TF activity scores in `grn_analysis/`
  - Motif enrichment analysis
  - Regulatory network edges (TF → target gene)
  - Comparative analysis (Control vs Mutant)
- **Dependencies**: 10_integration_prep.sh
- **Resources**: 16 CPUs, 64GB RAM, 12h
- **Key Feature**: Identifies disrupted regulatory networks

#### 12. **12_multiqc_report.sh** ✨ NEW
- **Purpose**: Generate comprehensive HTML quality report
- **Input**: All QC files and logs from previous steps
- **Output**: 
  - `multiqc_report/scATAC_QC_report.html`
  - Interactive quality dashboard
- **Dependencies**: 08_qc_metrics.sh (recommended)
- **Features**: 
  - Sample comparison tables
  - Fragment length distributions
  - TSS enrichment profiles
  - Processing pipeline overview

## Execution Instructions

### Run Core Pipeline (Required)
```bash
# Run remaining core scripts in order:
sbatch 03_validate_counts.sh
# Wait for completion, then:
sbatch 05_call_peaks.sh  
# Wait for completion, then:
sbatch 06_peak_cell_matrix.sh
```

### Run Advanced Analysis (Optional but Recommended)
```bash
# After step 6 completes, run these in parallel:
sbatch 07_create_bigwig.sh
sbatch 08_qc_metrics.sh
sbatch 09_dimensionality_reduction.sh
sbatch 10_integration_prep.sh

# After step 10 completes:
sbatch 11_gene_regulatory_analysis.sh

# Generate final report (run last):
sbatch 12_multiqc_report.sh
```

## Pipeline Dependencies Graph

```
FASTQ Files
     ↓
[01] Extract Barcodes
     ↓
[02] Test Barcodes
     ↓
┌────────────────┬────────────────┐
│ [03] Validate  │ [04] Chromap   │
│     Counts     │    Alignment   │
└────────────────┴────────┬───────┘
                          ↓
                 [05] Call Peaks
                          ↓
                [06] Peak Cell Matrix
                          ↓
    ┌─────────────────────┼─────────────────────┐
    ↓                     ↓                     ↓
[07] BigWig        [08] QC Metrics    [09] Dim Reduction
    ↓                     ↓                     ↓
    │              [10] Integration    [11] GRN Analysis
    │                     ↓                     ↓
    └─────────────→ [12] MultiQC Report ←──────┘
```

## Environment Requirements

- **Main environment**: `alignment_two` (conda environment)
- **BigWig environment**: `bigwig` (for script 07)
- **MACS2 environment**: `macs2_env` (for script 06)

## Key Output Directories

- `fragments/` - Fragment files from chromap alignment
- `peaks/` - Called peaks from MACS2
- `{sample}_peak_bc_matrix/` - 10X format count matrices
- `bigwig/` - Coverage tracks for genome browsers
- `qc_metrics/` - Quality control analysis results
- `dimensionality_reduction/` - LSI/PCA/UMAP results
- `integration_prep/` - Data prepared for scRNA integration
- `grn_analysis/` - Gene regulatory network analysis
- `multiqc_report/` - Final HTML quality report

## Sample Information

**Samples being processed:**
- `R26-Nestin-Ctrl-adult` (Control)
- `R26-Nestin-Mut-adult` (Mutant)

**Goal**: Identify disrupted gene regulatory networks in R26-Nestin mutant neural cells using integrated scATAC-seq and scRNA-seq analysis.

## Current Status

**Completed Steps:** ✅
- Barcode extraction and validation
- Chromap alignment and fragment generation

**Pending Core Steps:** ❌  
- Count validation (step 3)
- Peak calling (step 5)  
- Peak-cell matrix generation (step 6)

**Advanced Analysis:** 🚀
- All advanced scripts (7-12) are ready to run after core completion

## Notes

- Scripts 7-11 can run in parallel after step 6 completes
- Script 11 (GRN analysis) is specifically designed for your regulatory network analysis goals
- Script 12 generates the final comprehensive QC report
- All scripts are configured for SLURM job arrays processing both samples simultaneously