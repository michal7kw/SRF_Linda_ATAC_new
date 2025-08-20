#!/bin/bash

# Alternative tool installation script to handle conda dependency conflicts
# This script tries multiple approaches to install the required tools

set -e

echo "Installing bioinformatics tools (handling dependency conflicts)..."

# Method 1: Install tools one by one (avoiding MACS2 conflicts)
echo "Installing BWA and samtools first..."
conda install -c bioconda bwa samtools -y

# Method 2: Install MACS2 via pip (avoids conda conflicts)
echo "Installing MACS2 via pip..."
pip install MACS2

# Method 3: Install bowtie2 separately
echo "Installing bowtie2..."
conda install -c bioconda bowtie2 -y || echo "Bowtie2 installation optional - BWA will be used"

# Verify installations
echo ""
echo "Verifying installations..."

if command -v bwa &> /dev/null; then
    echo "✓ BWA installed: $(which bwa)"
    bwa 2>&1 | head -3
else
    echo "✗ BWA installation failed"
    exit 1
fi

if command -v samtools &> /dev/null; then
    echo "✓ samtools installed: $(which samtools)"
    samtools --version | head -1
else
    echo "✗ samtools installation failed"
    exit 1
fi

if command -v macs2 &> /dev/null; then
    echo "✓ MACS2 installed: $(which macs2)"
    macs2 --version
else
    echo "✗ MACS2 installation failed"
    exit 1
fi

if command -v bowtie2 &> /dev/null; then
    echo "✓ Bowtie2 installed: $(which bowtie2)"
    bowtie2 --version | head -1
else
    echo "! Bowtie2 not installed (optional - BWA will be used)"
fi

echo ""
echo "Tool installation completed successfully!"
echo "Core tools (BWA, samtools, MACS2) are ready for bulk ATAC processing."
echo ""
echo "Now you can run: sbatch run_bulk_atac_processing.sh"