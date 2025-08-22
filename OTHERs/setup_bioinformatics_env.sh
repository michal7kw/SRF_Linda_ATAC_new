#!/bin/bash

# Create dedicated bioinformatics conda environment with compatible Python version
# This avoids conflicts with your base Python 3.12 environment

set -e

echo "Setting up dedicated bioinformatics conda environment..."

# Initialize conda for this script
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh

# Environment name
ENV_NAME="atac_analysis"

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '$ENV_NAME' already exists. Removing it first..."
    conda env remove -n $ENV_NAME -y
fi

echo "Creating new conda environment '$ENV_NAME' with Python 3.10..."
conda create -n $ENV_NAME python=3.10 -y

echo "Activating environment and installing bioinformatics tools..."
# Use conda run to install packages in the new environment
conda run -n $ENV_NAME conda install -c bioconda -c conda-forge bwa samtools macs2 bowtie2 -y

echo "Verifying installations in the new environment..."
echo ""

# Test BWA
if conda run -n $ENV_NAME which bwa &> /dev/null; then
    echo "✓ BWA installed: $(conda run -n $ENV_NAME which bwa)"
    conda run -n $ENV_NAME bwa 2>&1 | head -3
else
    echo "✗ BWA installation failed"
    exit 1
fi

# Test samtools
if conda run -n $ENV_NAME which samtools &> /dev/null; then
    echo "✓ samtools installed: $(conda run -n $ENV_NAME which samtools)"
    conda run -n $ENV_NAME samtools --version | head -1
else
    echo "✗ samtools installation failed"
    exit 1
fi

# Test MACS2
if conda run -n $ENV_NAME which macs2 &> /dev/null; then
    echo "✓ MACS2 installed: $(conda run -n $ENV_NAME which macs2)"
    conda run -n $ENV_NAME macs2 --version
else
    echo "✗ MACS2 installation failed"
    exit 1
fi

# Test bowtie2
if conda run -n $ENV_NAME which bowtie2 &> /dev/null; then
    echo "✓ Bowtie2 installed: $(conda run -n $ENV_NAME which bowtie2)"
    conda run -n $ENV_NAME bowtie2 --version | head -1
else
    echo "! Bowtie2 installation failed (optional)"
fi

echo ""
echo "=========================================="
echo "✓ Bioinformatics environment setup complete!"
echo "Environment name: $ENV_NAME"
echo "Python version: $(conda run -n $ENV_NAME python --version)"
echo "=========================================="
echo ""
echo "To activate this environment manually:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To deactivate when done:"
echo "  conda deactivate"
echo ""
echo "The bulk processing script will automatically use this environment."
echo "You can now run: sbatch run_bulk_atac_processing.sh"