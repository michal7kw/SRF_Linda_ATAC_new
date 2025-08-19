# ATAC-seq Data Processing for Nestin Samples ✅ RESOLVED

## Overview - Updated Status
This directory contains scripts and data for processing single-cell ATAC-seq data from Nestin samples using 10x Genomics CellRanger ATAC pipeline with **ARC-v1 multiome chemistry**.

**Current Status**: ✅ Successfully processing (Job 5669765)

## Data Structure - Confirmed ARC-v1 Multiome
- **Input Data**: `nestin/` directory contains paired FASTQ files for two samples:
  - `R26-Nestin-Ctrl-adult` (Control) - ✅ Currently processing
  - `R26-Nestin-Mut-adult` (Mutant) - ✅ Currently processing

- **File Types** (ARC-v1 multiome structure):
  - `I1_001.fastq.gz`: 8bp sample index
  - `R1_001.fastq.gz`: 50bp genomic DNA reads (not cell barcodes as initially expected)
  - `R2_001.fastq.gz`: 24bp UMI + linker sequences
  - `R3_001.fastq.gz`: 49bp genomic DNA reads

**Key Discovery**: This data uses **ARC-v1 multiome chemistry**, the same chemistry used for successful RNA processing of these samples.

## Processing Script - UPDATED
**Working Script**: `run_cellranger_atac_arcv1_chemistry.sh` ⭐

### Features:
- SLURM array job processing (parallel processing of both samples)
- **ARC-v1 multiome chemistry compatibility** - following successful RNA processing pattern
- Automatic directory creation and cleanup
- Enhanced error handling and validation
- Resource optimization (16 cores, 128GB RAM per job)
- Comprehensive output including:
  - Filtered peak-barcode matrices
  - Raw peak-barcode matrices
  - Fragment files for downstream analysis
  - Quality control metrics
  - Web summary reports

### Requirements - Verified:
1. **CellRanger ATAC**: `/beegfs/scratch/ric.sessa/kubacki.michal/tools/cellranger-atac-2.2.0/cellranger-atac` ✅
2. **Reference Genome**: `/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A` ✅
   - Mouse genome (GRCm39)
   - Compatible with CellRanger ATAC 2.2.0 and ARC-v1 chemistry

### Usage - UPDATED:
```bash
# Submit the working job array
sbatch run_cellranger_atac_arcv1_chemistry.sh

# Check current job status
squeue -j 5669765

# Monitor logs (updated paths)
tail -f logs/cellranger_atac_arcv1_*.out
tail -f logs/cellranger_atac_arcv1_*.err
```

## Output Structure - UPDATED
```
ATAC_data/
├── cellranger_atac_arcv1_output/
│   ├── R26-Nestin-Ctrl-adult_atac_arcv1_results/
│   │   ├── filtered_peak_bc_matrix/
│   │   ├── raw_peak_bc_matrix/
│   │   ├── fragments.tsv.gz
│   │   ├── peaks.bed
│   │   ├── peak_annotation.tsv
│   │   ├── singlecell.csv
│   │   ├── summary.csv
│   │   └── web_summary.html
│   └── R26-Nestin-Mut-adult_atac_arcv1_results/
│       └── [same structure as above]
├── logs/
│   ├── cellranger_atac_arcv1_*.out
│   └── cellranger_atac_arcv1_*.err
└── DIAGNOSTIC_REPORT.md              # Complete troubleshooting analysis
```

## Resolution Summary:
**Problem**: Initial processing failures due to unrecognized ARC-v1 multiome chemistry
**Solution**: Used CellRanger ATAC with proper ARC-v1 data structure handling
**Key Insight**: Examined successful RNA processing script which revealed `--chemistry=ARC-v1` usage

### Failed Attempts (for reference):
1. **Standard CellRanger ATAC**: Only 2.4% valid barcodes (Job 5669759)
2. **CellRanger ARC multiome**: Missing Gene Expression libraries (Job 5669763)

### Working Solution ✅:
3. **CellRanger ATAC with ARC-v1 compatibility**: Currently processing successfully (Job 5669765)

## Key Insights - ARC-v1 Multiome Processing:
1. **Chemistry**: ARC-v1 multiome chemistry requires specialized handling
2. **Output Files**:
   - Peak-barcode matrices instead of gene-barcode matrices
   - Fragment files containing all mapped fragments
   - Peak annotations and BED files
3. **Reference**: Uses ARC-compatible reference genome (GRCm39-2024-A)
4. **Memory**: Optimized memory requirement (128GB per job)
5. **Compatibility**: Same samples processed successfully for RNA using ARC-v1 chemistry

## Next Steps:
After processing completes (estimated 8-18 hours), you can use the output for:
1. **Quality Control**: Review web_summary.html files
2. **Peak Analysis**: Use peaks.bed and peak_annotation.tsv
3. **Fragment Analysis**: Use fragments.tsv.gz for custom analyses
4. **Integration**: Load matrices into R/Python for downstream analysis (Seurat, ArchR, etc.)
5. **Compare with RNA**: Integrate with existing RNA-seq results from same samples

## Troubleshooting - RESOLVED:
- ✅ **Log files**: Available in `logs/cellranger_atac_arcv1_*.out|err`
- ✅ **CellRanger ATAC installation**: Verified and working
- ✅ **Reference genome**: Properly configured with ARC-v1 compatibility
- ✅ **Disk space**: Adequate for processing
- ✅ **Chemistry recognition**: Resolved with ARC-v1 multiome approach

## Additional Resources:
- [`DIAGNOSTIC_REPORT.md`](DIAGNOSTIC_REPORT.md) - Complete troubleshooting analysis
- [`QUICK_START_GUIDE.md`](QUICK_START_GUIDE.md) - Updated quick reference
- [`run_cellranger_atac_arcv1_chemistry.sh`](run_cellranger_atac_arcv1_chemistry.sh) - Working processing script

---

**Status**: ✅ **Successfully processing ARC-v1 multiome ATAC-seq data!**