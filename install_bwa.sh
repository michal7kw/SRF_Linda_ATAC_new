#!/bin/bash

# Quick BWA installation script for bulk ATAC processing
# This downloads and compiles BWA if it's not available

set -euo pipefail

echo "Installing BWA for bulk ATAC processing..."

# Installation directory
TOOLS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/tools"
BWA_DIR="$TOOLS_DIR/bwa-0.7.17"

# Create tools directory if it doesn't exist
mkdir -p "$TOOLS_DIR"
cd "$TOOLS_DIR"

# Check if BWA is already installed
if [[ -x "$BWA_DIR/bwa" ]]; then
    echo "BWA already installed at $BWA_DIR/bwa"
    exit 0
fi

# Download and compile BWA
echo "Downloading BWA 0.7.17..."
wget -O bwa-0.7.17.tar.bz2 "https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2"

echo "Extracting BWA..."
tar -xjf bwa-0.7.17.tar.bz2

echo "Compiling BWA..."
cd bwa-0.7.17
make

echo "BWA installation completed!"
echo "BWA executable location: $BWA_DIR/bwa"
echo ""
echo "To use BWA, add it to your PATH or update the script to use: $BWA_DIR/bwa"

# Clean up
cd "$TOOLS_DIR"
rm -f bwa-0.7.17.tar.bz2

echo "Installation complete!"