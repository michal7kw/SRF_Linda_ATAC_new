#!/bin/bash
#SBATCH --job-name=multimodal_coupling
#SBATCH --output=logs/multimodal_coupling_%j.out
#SBATCH --error=logs/multimodal_coupling_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=workq

# Multimodal scRNA-seq and scATAC-seq Coupling Analysis
# =====================================================

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Create or activate environment with required packages
if ! conda env list | grep -q "multimodal_analysis"; then
    echo "Creating multimodal_analysis conda environment..."
    conda create -n multimodal_analysis python=3.9 -y
    conda activate multimodal_analysis
    
    # Install core packages
    conda install -c conda-forge pandas numpy matplotlib seaborn scipy -y
    conda install -c conda-forge scikit-learn umap-learn -y
    pip install scanpy anndata h5py
    
    # Optional: Install additional packages for advanced analysis
    pip install pybedtools episcanpy
    
    echo "Environment created and packages installed"
else
    echo "Activating existing multimodal_analysis environment..."
    conda activate multimodal_analysis
fi

# Set up error handling
set -euo pipefail

# Create logs directory
mkdir -p logs
mkdir -p coupling_analysis_output
mkdir -p integration_analysis_output

echo "========================================="
echo "MULTIMODAL COUPLING ANALYSIS"
echo "Start time: $(date)"
echo "========================================="

# Step 1: Run initial coupling analysis
echo "Step 1: Running barcode coupling analysis..."
python multimodal_coupling_analysis.py > logs/coupling_analysis.log 2>&1

if [ $? -eq 0 ]; then
    echo "✓ Coupling analysis completed successfully"
else
    echo "✗ Coupling analysis failed - check logs/coupling_analysis.log"
    exit 1
fi

# Step 2: Run advanced integration if coupling was successful
if [ -d "coupling_analysis_output" ] && [ "$(ls -A coupling_analysis_output/*.h5ad 2>/dev/null)" ]; then
    echo "Step 2: Running advanced integration pipeline..."
    python advanced_integration_pipeline.py > logs/integration_analysis.log 2>&1
    
    if [ $? -eq 0 ]; then
        echo "✓ Integration analysis completed successfully"
    else
        echo "✗ Integration analysis failed - check logs/integration_analysis.log"
        echo "Coupling results are still available in coupling_analysis_output/"
    fi
else
    echo "⚠ No coupled datasets found - skipping integration step"
fi

echo "========================================="
echo "ANALYSIS SUMMARY"
echo "End time: $(date)"
echo "========================================="

# Generate summary report
echo "Creating analysis summary..."

SUMMARY_FILE="coupling_analysis_output/analysis_summary.txt"
{
    echo "Multimodal scRNA-seq/scATAC-seq Coupling Analysis Summary"
    echo "Generated: $(date)"
    echo "========================================================"
    echo ""
    
    echo "Input Data:"
    echo "- RNA data processed with CellRanger (ARC-v1 chemistry)"
    echo "- ATAC data processed with custom pipeline (chromap alignment)"
    echo ""
    
    if [ -f "coupling_analysis_output/overlap_analysis_results.csv" ]; then
        echo "Coupling Results:"
        echo "- Overlap analysis completed"
        echo "- Results saved in: coupling_analysis_output/"
        
        # Count coupled datasets
        COUPLED_COUNT=$(ls coupling_analysis_output/*_coupled.h5ad 2>/dev/null | wc -l)
        echo "- Coupled datasets created: $COUPLED_COUNT"
        echo ""
        
        if [ -f "coupling_analysis_output/overlap_analysis_results.csv" ]; then
            echo "Barcode Overlap Summary:"
            head -1 "coupling_analysis_output/overlap_analysis_results.csv"
            tail -n +2 "coupling_analysis_output/overlap_analysis_results.csv" | while read line; do
                echo "  $line"
            done
            echo ""
        fi
    fi
    
    if [ -d "integration_analysis_output" ]; then
        echo "Integration Results:"
        INTEGRATED_COUNT=$(ls integration_analysis_output/*_integrated.h5ad 2>/dev/null | wc -l)
        echo "- Integrated datasets: $INTEGRATED_COUNT"
        echo "- Plots saved in: integration_analysis_output/"
        echo ""
    fi
    
    echo "Next Steps:"
    echo "1. Review QC plots to assess coupling quality"
    echo "2. Load .h5ad files in Python/R for downstream analysis"
    echo "3. Run differential analysis between conditions"
    echo "4. Perform trajectory analysis or other specialized analyses"
    echo ""
    
    echo "Files Generated:"
    echo "- coupling_analysis_output/: Basic coupling results"
    echo "- integration_analysis_output/: Advanced integration results" 
    echo "- logs/: Analysis logs and debugging information"
    echo "- seurat_multimodal_analysis.R: Template for Seurat analysis"
    
} > "$SUMMARY_FILE"

echo "Analysis summary saved to: $SUMMARY_FILE"
echo ""

# Display key results
if [ -f "$SUMMARY_FILE" ]; then
    echo "=== ANALYSIS SUMMARY ==="
    cat "$SUMMARY_FILE"
fi

echo ""
echo "Analysis complete! Check the following directories:"
echo "- coupling_analysis_output/ : Coupled datasets and QC plots"
echo "- integration_analysis_output/ : Integrated analysis results"
echo "- logs/ : Detailed analysis logs"