#!/bin/bash
#SBATCH --job-name=validate_counts
#SBATCH --output=logs/validate_counts_%A_%a.out
#SBATCH --error=logs/validate_counts_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 3: Validating read counts for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "ERROR: Extracted barcode file not found: $EXTRACTED_BC_FILE"
    echo "Please run step 1 (01_extract_barcodes.sh) first"
    exit 1
fi

echo "DEBUG: Counting reads in each file..."
R1_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" | wc -l | awk '{print $1/4}')
R3_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" | wc -l | awk '{print $1/4}') 
BC_COUNT=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')

echo "Read counts:"
echo "  R1 (genomic, forward): $R1_COUNT"
echo "  R3 (genomic, reverse): $R3_COUNT" 
echo "  R2 (extracted barcodes): $BC_COUNT"

# Check if all counts match
if [[ "$R1_COUNT" -ne "$R3_COUNT" ]]; then
    echo "ERROR: R1 and R3 read counts don't match ($R1_COUNT vs $R3_COUNT)"
    exit 1
fi

if [[ "$R1_COUNT" -ne "$BC_COUNT" ]]; then
    echo "ERROR: R1 and barcode counts don't match ($R1_COUNT vs $BC_COUNT)"
    echo "This suggests the R2 barcode file was truncated or has different structure"
    echo "Attempting to fix by re-extracting barcodes..."
    
    # Re-extract barcodes ensuring same read count as R1
    echo "DEBUG: Re-extracting barcodes with read count validation..."
    ORIG_R2_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | wc -l | awk '{print $1/4}')
    echo "DEBUG: Original R2 read count: $ORIG_R2_COUNT"
    
    if [[ "$ORIG_R2_COUNT" -ne "$R1_COUNT" ]]; then
        echo "ERROR: Original R2 file has different read count than R1 ($ORIG_R2_COUNT vs $R1_COUNT)"
        echo "This indicates corrupted or incomplete FASTQ files"
        exit 1
    fi
    
    # Re-extract with validation
    rm -f "$EXTRACTED_BC_FILE"
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
        awk '{
            if(NR%4==1) print $0;
            else if(NR%4==2) print substr($0,9,16);
            else if(NR%4==3) print $0;
            else if(NR%4==0) print substr($0,9,16);
        }' | gzip > "$EXTRACTED_BC_FILE"
    
    # Recount after extraction
    BC_COUNT_NEW=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')
    echo "DEBUG: New barcode count after re-extraction: $BC_COUNT_NEW"
    
    if [[ "$BC_COUNT_NEW" -ne "$R1_COUNT" ]]; then
        echo "ERROR: Re-extraction failed, counts still don't match ($BC_COUNT_NEW vs $R1_COUNT)"
        exit 1
    fi
    
    echo "DEBUG: Re-extraction successful, all counts now match"
    BC_COUNT=$BC_COUNT_NEW
fi

echo "✓ All file read counts match: $R1_COUNT reads per file"

# Save validation results
mkdir -p "$OUTPUT_DIR/qc"
cat > "$OUTPUT_DIR/qc/${SAMPLE}_read_count_validation.txt" << EOF
SAMPLE=$SAMPLE
R1_COUNT=$R1_COUNT
R3_COUNT=$R3_COUNT
BC_COUNT=$BC_COUNT
VALIDATION_STATUS=PASS
VALIDATED_AT=$(date)
EOF

echo "Validation results saved to: $OUTPUT_DIR/qc/${SAMPLE}_read_count_validation.txt"

echo "========================================="
echo "Step 3 complete for $SAMPLE"
echo "All read counts validated: $R1_COUNT reads"
echo "End time: $(date)"
echo "========================================="