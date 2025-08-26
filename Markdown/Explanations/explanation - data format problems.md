---
created: 2025-08-26T08:16
updated: 2025-08-26T08:17
---
# Detailed Explanation of the ATAC-seq Data Processing Problem

This document provides a detailed explanation of the issues encountered while processing the ATAC-seq data for the Nestin samples (R26-Nestin-Ctrl-adult and R26-Nestin-Mut-adult).

## The Core Problem: Mismatch Between Data and Expected Format

The fundamental issue is a mismatch between the provided FASTQ data and the expected input format for the `cellranger-arc` pipeline, which is the standard tool for processing 10x Genomics Multiome ATAC + Gene Expression data. This mismatch manifests in two key areas:

1.  **Incorrect Read Structure:** The FASTQ files do not follow the standard read structure for ARC-v1 chemistry.
2.  **Low Barcode Recovery:** As a result of the incorrect read structure, `cellranger-arc` is unable to identify a sufficient number of valid 10x barcodes, leading to pipeline failure.

### Expected vs. Actual Read Structure

The following table details the differences between the expected and actual read structures:

| Read | Expected 10x Multiome ATAC + GEX (ARC-v1) | Actual Observed Data |
| :--- | :--- | :--- |
| **I1** | 8bp (Sample Index) | 8bp (Sample Index) |
| **R1** | 16bp (10x Barcode) | 50bp (Appears to be Read 1 - Genomic DNA) |
| **R2** | 50bp (Read 1 - Genomic DNA) | 24bp (Potentially contains the 10x Barcode and UMI) |
| **R3** | 49bp (Read 2 - Genomic DNA) | 49bp (Appears to be Read 2 - Genomic DNA) |

As the table illustrates, the 10x Barcode is expected to be in `R1`, but in the provided data, `R1` appears to contain genomic DNA. The 10x Barcode is likely located within the 24bp `R2` read, but this is not the standard location for `cellranger-arc`.

### Inconsistency with Documentation

The `cellranger-arc` documentation, in the "Specifying Input FASTQ Files for cellranger-arc count" section, describes the expected ATAC FASTQ format as:

*   `I1`: Dual index i7 read (optional)
*   `R1`: Read 1
*   `R2`: Dual index i5 read
*   `R3`: Read 2

The observed data does not match this format. The `R2` read in the provided data is 24bp and likely contains the barcode, while the documentation states that `R2` should be the dual index i5 read. This further confirms that the data is in a non-standard format that is incompatible with the `cellranger-arc` pipeline.

### Impact on Data Processing

This discrepancy in read structure has the following consequences:

*   **`cellranger-arc` cannot find the barcodes:** The pipeline is hard-coded to look for the 16bp 10x Barcode in the `R1` read. Since the barcodes are not in the expected location, the software fails to identify them.
*   **Low Barcode Recovery Rate:** The error message "1.4% (< 10%) of read pairs have a valid 10x barcode" is a direct result of this issue. The pipeline is unable to find the vast majority of barcodes, leading to a recovery rate far below the acceptable threshold.
*   **Pipeline Failure:** Because the barcode recovery rate is so low, `cellranger-arc` concludes that there is a problem with the data (e.g., poor sequencing quality, sample mixup, or incorrect pipeline) and terminates the analysis.

## Conclusion

The root cause of the data processing failure is the non-standard read structure of the provided FASTQ files. To resolve this issue, we need clarification from the sequencing facility on the exact library preparation and sequencing protocol used. This information will allow us to either adapt our processing pipeline to accommodate the non-standard read structure or to request that the data be re-sequenced in the correct format.