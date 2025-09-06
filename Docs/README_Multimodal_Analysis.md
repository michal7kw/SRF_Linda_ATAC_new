# Multimodal scRNA-seq and scATAC-seq Coupling Analysis

This repository contains a comprehensive pipeline for coupling and integrating single-cell RNA-seq and ATAC-seq data from multimodal experiments where the same nuclei were sequenced in separate runs.

## Overview

Your multimodal experiment setup:
- **RNA data**: Processed with CellRanger using ARC-v1 chemistry
- **ATAC data**: Custom processing pipeline with chromap alignment
- **Challenge**: Different barcode processing approaches need reconciliation
- **Goal**: Couple data from the same nuclei for integrated analysis

## Data Structure Analysis

### RNA Processing (CellRanger)
- **Script**: `run_cellranger_cluster_only_counts.sh`
- **Chemistry**: ARC-v1 (multiome compatible)
- **Samples**: R26_Emx1_Ctrl_adult, R26_Emx1_Mut_adult, R26_Nestin_Ctrl_adult, R26_Nestin_Mut_adult
- **Output**: Standard 10X matrices with 16bp barcodes + `-1` suffix

### ATAC Processing (Custom Pipeline)
- **Scripts**: `01_extract_barcodes.sh` → `02_test_barcodes.sh` → `03_validate_counts.sh` → `04_chromap_alignment.sh`
- **Barcode extraction**: Rightmost 16bp from 24bp R2 reads (positions 9-24)
- **Samples**: R26-Nestin-Ctrl-adult, R26-Nestin-Mut-adult
- **Output**: Chromap fragment files with 16bp barcodes (no suffix)

## Pipeline Components

### 1. Multimodal Coupling Analysis (`multimodal_coupling_analysis.py`)

**Main class**: `MultimodalCoupler`

**Key functions**:
- `find_rna_data_paths()`: Locate CellRanger output matrices
- `find_atac_data_paths()`: Locate ATAC fragment files  
- `load_rna_data()`: Load 10X RNA matrices with scanpy
- `load_atac_fragments()`: Parse ATAC fragment files and extract barcode statistics
- `analyze_barcode_overlap()`: Calculate overlap between RNA and ATAC barcodes
- `create_coupled_dataset()`: Generate datasets containing only overlapping cells
- `run_coupling_analysis()`: Execute complete coupling pipeline

**Output**:
- Coupled AnnData objects (`.h5ad` files)
- Barcode overlap statistics
- QC plots showing coupling quality

### 2. Advanced Integration Pipeline (`advanced_integration_pipeline.py`)

**Main class**: `AdvancedMultimodalIntegrator`

**Key functions**:
- `call_peaks_from_fragments()`: Simple peak calling from ATAC fragments
- `create_peak_cell_matrix()`: Generate peak-by-cell count matrix
- `calculate_gene_activity_scores()`: Compute gene activity from peak accessibility
- `integrate_multimodal_data()`: Combine RNA and ATAC modalities
- `run_integration_analysis()`: Dimensionality reduction and clustering
- `generate_integration_plots()`: Comprehensive visualization

**Methods supported**:
- Concatenation approach (RNA + gene activity scores)
- WNN-style integration (simplified implementation)
- Extensible framework for additional methods

## Usage

### Quick Start

1. **Run the complete pipeline**:
```bash
# Make script executable
chmod +x run_multimodal_analysis.sh

# Submit to cluster
sbatch run_multimodal_analysis.sh
```

2. **Monitor progress**:
```bash
# Check job status
squeue -u $USER

# Monitor logs
tail -f logs/multimodal_coupling_*.out
```

### Manual Execution

1. **Coupling analysis only**:
```bash
# Activate environment
conda activate multimodal_analysis  # or create new environment

# Run coupling
python multimodal_coupling_analysis.py
```

2. **Advanced integration** (after coupling):
```bash
python advanced_integration_pipeline.py
```

### Python Environment Setup

```bash
# Create environment
conda create -n multimodal_analysis python=3.9 -y
conda activate multimodal_analysis

# Install packages
conda install -c conda-forge pandas numpy matplotlib seaborn scipy -y
conda install -c conda-forge scikit-learn umap-learn -y
pip install scanpy anndata h5py pybedtools episcanpy
```

## Key Outputs

### Coupling Analysis Results
- **`coupling_analysis_output/`**:
  - `*_coupled.h5ad`: Coupled datasets with both RNA and ATAC info
  - `*_barcodes.txt`: Lists of overlapping cell barcodes  
  - `overlap_analysis_results.csv`: Statistical summary of barcode overlaps
  - `barcode_overlap_summary.png`: Visual summary of coupling quality
  - `qc_*.png`: Per-sample quality control plots

### Integration Analysis Results
- **`integration_analysis_output/`**:
  - `*_integrated.h5ad`: Fully integrated multimodal datasets
  - `*_peaks.csv`: Called ATAC peaks
  - `*_gene_activity.h5ad`: Gene activity score matrices
  - `integration_*.png`: UMAP plots and clustering results
  - `seurat_multimodal_analysis.R`: Template R script for Seurat analysis

## Data Flow

```
RNA Data (10X matrices) ─┐
                         ├─→ Barcode Matching ─→ Coupled Datasets ─→ Integration ─→ Analysis
ATAC Data (fragments) ───┘
```

1. **Load data**: RNA matrices + ATAC fragments
2. **Extract barcodes**: Clean barcode sequences from both modalities
3. **Find overlaps**: Identify cells present in both datasets
4. **Create coupled data**: Filter to overlapping cells only
5. **Call peaks**: Generate ATAC peak-cell matrices
6. **Gene activity**: Convert peaks to gene-level scores
7. **Integration**: Combine modalities for joint analysis
8. **Analysis**: Clustering, UMAP, differential analysis

## Sample Matching Logic

The pipeline automatically matches samples based on naming patterns:

| RNA Sample | ATAC Sample | Match Status |
|------------|-------------|--------------|
| R26_Nestin_Ctrl_adult | R26-Nestin-Ctrl-adult | ✓ Matched |
| R26_Nestin_Mut_adult | R26-Nestin-Mut-adult | ✓ Matched |
| R26_Emx1_Ctrl_adult | R26-Emx1-Ctrl-adult | ⚠ Check data availability |
| R26_Emx1_Mut_adult | R26-Emx1-Mut-adult | ⚠ Check data availability |

## Troubleshooting

### Common Issues

1. **Low barcode overlap**:
   - Check barcode extraction parameters in ATAC processing
   - Verify whitelist compatibility between RNA and ATAC
   - Consider barcode error correction

2. **Memory issues**:
   - Increase SLURM memory allocation
   - Use fragment file subsampling for initial analysis
   - Process samples individually

3. **Missing data**:
   - Verify file paths in base directory structure
   - Check CellRanger output completeness
   - Ensure ATAC fragment files are properly indexed

### Expected Results

- **Good coupling**: 20-80% barcode overlap between modalities
- **Typical cell numbers**: 1,000-10,000 coupled cells per sample
- **Quality metrics**: >500 genes/cell (RNA), >1,000 fragments/cell (ATAC)

## Next Steps for Analysis

1. **Quality control**: Review coupling efficiency and cell quality metrics
2. **Differential analysis**: Compare gene expression and chromatin accessibility between conditions
3. **Trajectory analysis**: Pseudotime analysis using integrated data
4. **Regulatory analysis**: Link accessible chromatin regions to gene expression
5. **Cell type identification**: Cluster analysis and marker gene identification

## Integration with Other Tools

### Seurat (R)
```r
# Load integrated data
library(Seurat)
library(Signac)

# Convert from Python/scanpy format if needed
# Use sceasy or similar conversion tools
```

### ArchR (R)
```r
# For advanced ATAC-seq analysis
library(ArchR)

# Import peak matrices and fragment files
# Run peak calling, gene activity scoring, footprinting
```

### SCENIC+ (Python)
```python
# For gene regulatory network analysis
import pyscenic

# Use integrated data for enhancer-gene linking
# Infer transcription factor activities
```

## Contact and Support

For questions about this pipeline or multimodal analysis approaches:

1. Check logs in `logs/` directory for detailed error messages
2. Review QC plots for data quality assessment  
3. Consult scanpy/Seurat documentation for advanced analysis options

## References

Key methods and tools used:
- **Scanpy**: Wolf et al., Genome Biology 2018
- **10X Genomics multiome**: Single Cell Gene Expression with Feature Barcoding
- **Chromap**: Zhang et al., Nature Communications 2021
- **WNN integration**: Hao et al., Cell 2021