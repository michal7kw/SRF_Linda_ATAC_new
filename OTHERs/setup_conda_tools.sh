#!/bin/bash

# Setup script to install required bioinformatics tools via conda
# Run this before using the bulk ATAC processing script

echo "Installing bioinformatics tools via conda..."

# Install BWA for alignment
echo "Installing BWA..."
conda install -c bioconda bwa -y

# Install other useful tools for ATAC-seq analysis
echo "Installing additional ATAC-seq tools..."
conda install -c bioconda samtools -y
conda install -c bioconda macs2 -y
conda install -c bioconda bowtie2 -y  # backup aligner

# Verify installations
echo ""
echo "Verifying installations..."

if command -v bwa &> /dev/null; then
    echo "✓ BWA installed: $(which bwa)"
    bwa 2>&1 | head -3
else
    echo "✗ BWA installation failed"
fi

if command -v samtools &> /dev/null; then
    echo "✓ samtools installed: $(which samtools)"
    samtools --version | head -1
else
    echo "✗ samtools installation failed"
fi

if command -v macs2 &> /dev/null; then
    echo "✓ MACS2 installed: $(which macs2)"
    macs2 --version
else
    echo "✗ MACS2 installation failed"
fi

if command -v bowtie2 &> /dev/null; then
    echo "✓ Bowtie2 installed: $(which bowtie2)"
    bowtie2 --version | head -1
else
    echo "✗ Bowtie2 installation failed"
fi

echo ""
echo "Setup complete! You can now run the bulk ATAC processing script."
echo "Usage: sbatch run_bulk_atac_processing.sh"