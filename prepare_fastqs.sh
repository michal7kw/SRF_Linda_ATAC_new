#!/bin/bash
# This script prepares the FASTQ files for CellRanger ARC by creating symbolic links
# with the correct naming convention (SampleName_S1_L001_R1_001.fastq.gz).

set -euo pipefail

# Configuration
SOURCE_DIR="ATAC_data/nestin"
TARGET_DIR="ATAC_data/nestin_prepared"
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# Create the target directory
mkdir -p "$TARGET_DIR"
echo "Created target directory: $TARGET_DIR"

# Create symbolic links for each sample and read type
for SAMPLE in "${SAMPLES[@]}"; do
    for READ in "I1" "R1" "R2" "R3"; do
        SOURCE_FILE="$SOURCE_DIR/${SAMPLE}_${READ}_001.fastq.gz"
        TARGET_FILE="$TARGET_DIR/${SAMPLE}_S1_L001_${READ}_001.fastq.gz"
        
        # Get the absolute path of the source file for the symbolic link
        ABS_SOURCE_FILE="$(cd "$SOURCE_DIR" && realpath "${SAMPLE}_${READ}_001.fastq.gz")"
        
        if [ -f "$ABS_SOURCE_FILE" ]; then
            ln -sf "$ABS_SOURCE_FILE" "$TARGET_FILE"
            echo "Linked: $TARGET_FILE -> $ABS_SOURCE_FILE"
        else
            echo "WARNING: Source file not found: $ABS_SOURCE_FILE"
        fi
    done
done

echo "FASTQ file preparation complete."