---
created: 2025-08-26T08:16
updated: 2025-08-26T08:17
---
# Explanation of the 10x Barcode and UMI Extraction Step

This document explains the necessity of extracting the 10x Barcode and UMI from the `R2` read and how to perform this extraction.

## The Problem: Non-Standard Read Structure

As detailed in the `problem_explanation.md` document, the provided FASTQ files have a non-standard read structure. The 10x Barcode, which is expected to be in the `R1` read, is instead embedded within the `R2` read.

**Observed R2 Read Structure:**

```
[CAGACG][16bp barcode][XX]
  └─6bp─┘ └────16bp────┘└2bp┘
  Total: 24bp
```

The 16bp barcode is flanked by a 6bp linker sequence (`CAGACG`) and 2bp of genomic DNA. This non-standard structure prevents `cellranger-arc` from identifying the barcodes, leading to pipeline failure.

## The Solution: Barcode Extraction

To resolve this issue, we must extract the 16bp barcode from the `R2` read and create a new FASTQ file that contains only the barcode. This new FASTQ file will then be used as the `R1` read for `cellranger-arc`.

### Extraction Method

The following `awk` command can be used to extract the 16bp barcode from the `R2` read:

```bash
awk '{if(NR%4==2) print substr($0,7,16)}'
```

**Explanation:**

*   `awk`: A powerful command-line tool for text processing.
*   `'{if(NR%4==2) ... }'`: This is an `awk` script that is executed for each line of the input file.
    *   `NR`: The current record (line) number.
    *   `NR%4==2`: This condition is true for every second line of a 4-line FASTQ record (the sequence line).
*   `print substr($0,7,16)`: If the condition is true, this command is executed.
    *   `print`: Prints the output to the console.
    *   `substr($0,7,16)`: Extracts a substring from the current line (`$0`).
        *   `7`: The starting position of the substring (the 7th character).
        *   `16`: The length of the substring (16 characters).

This command will effectively extract the 16bp barcode from each `R2` read and print it to a new file.

### Implementation in the Script

The `run_cellranger_atac_arcv1_chemistry_corrected.sh` script has been modified to incorporate this barcode extraction step. The script now performs the following actions:

1.  **Creates a new `PROCESSED_DATA_DIR`:** This directory will store the extracted barcode FASTQ files.
2.  **Extracts the barcode:** The `awk` command is used to extract the 16bp barcode from the `R2` read and create a new FASTQ file named `[SAMPLE_NAME]_R2_16bp.fastq.gz`.
3.  **Renames the files:** The script renames the FASTQ files to match the expected format for `cellranger-arc`. The extracted barcode FASTQ file is used as the `R1` read.

By performing this barcode extraction and file renaming, we can provide `cellranger-arc` with the correctly formatted input files, allowing the pipeline to proceed with the analysis.