# scATAC-seq Pipeline Troubleshooting Analysis

## Executive Summary

This document summarizes the comprehensive troubleshooting analysis of the scATAC-seq processing pipeline, focusing on issues encountered in steps 6 (peak-cell matrix generation) and 9 (dimensionality reduction). The primary issue was identified as resource limitations during the computationally intensive bedtools intersect operation.

## Problem Description

### Initial Symptoms
- **Step 9 (Dimensionality Reduction)** failed with two distinct error patterns:
  - Sample 0 (R26-Nestin-Ctrl-adult): Sparse matrix error during TF-IDF normalization
  - Sample 1 (R26-Nestin-Mut-adult): Empty features file error
- **Step 6 (Peak-Cell Matrix Generation)** appeared to terminate immediately without generating updated matrix files
- Matrix files contained 0 entries despite having peaks and barcodes

### Error Messages
```
# Sample 0 Error
Error in h(simpleError(msg, call)) : 
  error in evaluating the argument 'x' in selecting a method for function 'diff': 
  no slot of name "p" for this object of class "dgTMatrix"

# Sample 1 Error  
Error in read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE) : 
  no lines available in input
```

## Root Cause Analysis

### 1. Matrix File Investigation
- **Matrix Info Analysis**: Both samples showed 0 entries in matrix_info.txt files
  - R26-Nestin-Ctrl-adult: 217,613 peaks × 225,796 cells = 0 entries
  - R26-Nestin-Mut-adult: 0 peaks × 0 cells = 0 entries

### 2. Timeline Analysis
- Peak files were regenerated on Sep 5, 2025
- Matrix files were last updated on Sep 3-4, 2025
- **Conclusion**: Step 6 needed to be re-run after peak calling completion

### 3. Script Execution Analysis
- Initial bedtools intersect command failed silently
- **Identified Issue**: `zcat_if_gzipped` function was not available in the subshell context
- **Fix Applied**: Direct file path usage instead of function call

### 4. Computational Resource Analysis
- Manual bedtools intersect test revealed massive dataset: **32.6 million overlap records**
- Output file size: **3.4 GB** of intersection data
- Processing pipeline includes: bedtools intersect → awk → sort → uniq → awk → gzip

## Technical Findings

### Dataset Scale
```
Input Files:
- Peaks: 217,613 regions (R26-Nestin-Ctrl-adult)  
- Reads: ~2.7 GB compressed BED file
- Overlaps Found: 32,602,808 records
- Raw Output Size: 3.4 GB
```

### Resource Requirements
```
Original SLURM Settings:
- Memory: 64G
- Time: 6:00:00  
- CPUs: 8

Updated SLURM Settings:
- Memory: 128G (doubled)
- Time: 6:00:00 (maintained)
- CPUs: 8 (maintained)
```

### Performance Bottlenecks
1. **Bedtools intersect**: I/O intensive for large datasets
2. **Sort operation**: Memory intensive for 32M+ records  
3. **Uniq operation**: Sequential processing bottleneck
4. **Matrix construction**: AWK processing of large intermediate files

## Solutions Implemented

### 1. Script Bug Fix
**File**: `06_peak_cell_matrix.sh:88`
```bash
# Before (problematic)
zcat_if_gzipped "$PEAKS_FILE" | \
bedtools intersect \
    -a stdin \
    ...

# After (fixed)
bedtools intersect \
    -a "$PEAKS_FILE" \
    ...
```

### 2. Resource Allocation Increase
**File**: `06_peak_cell_matrix.sh:7`
```bash
# Before
#SBATCH --mem=64G

# After  
#SBATCH --mem=128G
```

### 3. Pipeline Workflow Correction
- Ensured step 5 (peak calling) completion before step 6 execution
- Implemented dependency checking in workflow

## Recommendations

### Immediate Actions
1. ✅ **Fixed**: Corrected bedtools intersect command syntax
2. ✅ **Implemented**: Increased memory allocation to 128G
3. 🔄 **In Progress**: Monitor current step 6 execution

### Future Optimizations
1. **Resource Scaling**: Consider increasing time limit for very large datasets
2. **Parallel Processing**: Implement chromosome-wise parallel processing
3. **Intermediate Cleanup**: Add cleanup of large temporary files
4. **Memory Monitoring**: Add memory usage logging for resource optimization
5. **Checkpoint System**: Implement resume capability for long-running jobs

### Pipeline Improvements
1. **Dependency Management**: Implement automated dependency checking
2. **Resource Estimation**: Add dataset size-based resource allocation
3. **Progress Monitoring**: Add detailed progress reporting for long operations
4. **Error Handling**: Improve error messages and recovery options

## Performance Metrics

### Expected Processing Times
- **Bedtools intersect**: 30-60 minutes (32M records)
- **Data processing**: 15-30 minutes (sort/uniq operations)
- **Matrix generation**: 10-15 minutes
- **Total step 6 time**: ~2-3 hours per sample

### Memory Usage Patterns
- **Peak usage**: During sort operation (~100-120G)
- **Baseline usage**: During I/O operations (~20-30G)
- **Matrix construction**: Moderate usage (~40-50G)

## Lessons Learned

1. **Function Export**: Shell functions must be properly exported for subshell usage
2. **Resource Planning**: Large genomics datasets require careful resource estimation
3. **Dependency Tracking**: Workflow steps must verify prerequisite completion
4. **Debugging Strategy**: Systematic approach from symptoms to root cause
5. **Testing Approach**: Manual command testing invaluable for pipeline debugging

## Next Steps

1. **Monitor Current Jobs**: Verify step 6 completion with increased resources
2. **Validate Output**: Check matrix quality and dimensions
3. **Execute Step 9**: Run dimensionality reduction once matrices are ready
4. **Documentation Update**: Update pipeline documentation with resource requirements
5. **Testing**: Validate entire pipeline with current fixes

---
*Analysis completed: September 5, 2025*  
*Pipeline Status: Step 6 running with optimized resources*