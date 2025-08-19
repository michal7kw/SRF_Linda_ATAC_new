#!/bin/bash
# Install FastQC

TOOLS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/tools"
mkdir -p $TOOLS_DIR
cd $TOOLS_DIR

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod +x FastQC/fastqc