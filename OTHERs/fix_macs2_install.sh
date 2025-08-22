#!/bin/bash

# This script fixes the MACS2 installation in the atac_analysis conda environment
# by reinstalling it from source to resolve library conflicts.

set -e

ENV_NAME="atac_analysis"

echo "Fixing MACS2 installation in conda environment: $ENV_NAME"

# Initialize conda for this script
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh

# Activate the environment
conda activate $ENV_NAME

echo "Environment activated. Python version: $(python --version)"

# Uninstall the potentially broken MACS2 package
echo "Uninstalling existing MACS2..."
pip uninstall -y MACS2

# Reinstall MACS2 using pip with --no-binary flag to force compilation
echo "Reinstalling MACS2 from source using pip..."
pip install --no-binary :all: MACS2

echo "Verifying MACS2 installation..."

if command -v macs2 &> /dev/null; then
    echo "✓ MACS2 reinstalled successfully: $(which macs2)"
    macs2 --version
else
    echo "✗ MACS2 reinstallation failed."
    exit 1
fi

echo ""
echo "======================================================"
echo "MACS2 has been fixed in the '$ENV_NAME' environment."
echo "You can now re-run the failed peak calling step or the full script."
echo ""
echo "To re-run the full script for the emx1 sample:"
echo "  sbatch ATAC_data/run_bulk_atac_emx1.sh"
echo "======================================================"

conda deactivate