#!/usr/bin/env python3
"""
Multimodal scRNA-seq and scATAC-seq Data Coupling Analysis
=========================================================

This script provides comprehensive tools for coupling single-cell RNA and ATAC-seq data
from multimodal experiments where the same nuclei were sequenced across different runs.

Key functionalities:
1. Load and preprocess both RNA and ATAC data
2. Analyze barcode overlap and quality
3. Create coupled datasets for downstream analysis
4. Generate comprehensive QC reports

Author: Claude Code Assistant
Date: September 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
from scipy import sparse
import h5py
import gzip
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

class MultimodalCoupler:
    """
    Main class for coupling scRNA and scATAC data from multimodal experiments
    """
    
    def __init__(self, base_dir="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda"):
        self.base_dir = Path(base_dir)
        self.rna_data = {}
        self.atac_data = {}
        self.coupled_data = {}
        self.results = {}
        
        # Define sample mapping between RNA and ATAC naming conventions
        self.sample_mapping = {
            'RNA': {
                'R26_Nestin_Ctrl_adult': ['R26-Nestin-Ctrl-adult', 'R26_Nestin_Ctrl_adult'],
                'R26_Nestin_Mut_adult': ['R26-Nestin-Mut-adult', 'R26_Nestin_Mut_adult'],
                'R26_Emx1_Ctrl_adult': ['R26-Emx1-Ctrl-adult', 'R26_Emx1_Ctrl_adult'], 
                'R26_Emx1_Mut_adult': ['R26-Emx1-Mut-adult', 'R26_Emx1_Mut_adult']
            }
        }
        
        print(f"Initialized MultimodalCoupler with base directory: {self.base_dir}")
    
    def find_rna_data_paths(self):
        """Find RNA data matrices from CellRanger output"""
        rna_paths = {}
        
        # Search in multiple potential locations
        search_dirs = [
            self.base_dir / "SRF_Linda_RNA" / "cellranger_final_count_data_trans",
            self.base_dir / "SRF_Spatial_segmentation" / "DATA" / "RNA_counts"
        ]
        
        for search_dir in search_dirs:
            if search_dir.exists():
                print(f"Searching for RNA data in: {search_dir}")
                
                for sample_dir in search_dir.iterdir():
                    if sample_dir.is_dir():
                        # Look for filtered feature matrix
                        matrix_path = sample_dir / "outs" / "filtered_feature_bc_matrix"
                        if matrix_path.exists():
                            sample_name = sample_dir.name.replace("cellranger_counts_", "")
                            rna_paths[sample_name] = matrix_path
                            print(f"Found RNA data for {sample_name}: {matrix_path}")
        
        return rna_paths
    
    def find_atac_data_paths(self):
        """Find ATAC fragment files"""
        atac_paths = {}
        
        atac_dir = self.base_dir / "ATAC_data" / "chromap_final_output" / "fragments"
        
        if atac_dir.exists():
            print(f"Searching for ATAC data in: {atac_dir}")
            
            for fragment_file in atac_dir.glob("*_fragments.tsv.gz"):
                sample_name = fragment_file.stem.replace("_fragments.tsv", "")
                atac_paths[sample_name] = fragment_file
                print(f"Found ATAC data for {sample_name}: {fragment_file}")
        
        return atac_paths
    
    def load_rna_data(self, sample_name, matrix_path):
        """Load 10X RNA data from CellRanger output"""
        print(f"Loading RNA data for {sample_name}...")
        
        try:
            # Load using scanpy
            adata = sc.read_10x_mtx(
                matrix_path,
                var_names='gene_symbols',
                cache=True,
                gex_only=True
            )
            
            # Make variable names unique
            adata.var_names_unique()
            
            # Add sample information
            adata.obs['sample'] = sample_name
            adata.obs['modality'] = 'RNA'
            
            # Extract clean barcodes (remove -1 suffix if present)
            adata.obs['barcode_clean'] = adata.obs.index.str.split('-').str[0]
            
            print(f"  Loaded {adata.n_obs} cells x {adata.n_vars} genes")
            print(f"  Sample barcodes: {adata.obs['barcode_clean'].head().tolist()}")
            
            return adata
            
        except Exception as e:
            print(f"Error loading RNA data for {sample_name}: {e}")
            return None
    
    def load_atac_fragments(self, sample_name, fragment_path):
        """Load ATAC fragment file and extract barcode information"""
        print(f"Loading ATAC fragments for {sample_name}...")
        
        try:
            # Read fragment file (compressed)
            fragments = []
            barcode_counts = {}
            
            with gzip.open(fragment_path, 'rt') as f:
                for i, line in enumerate(f):
                    if i % 100000 == 0 and i > 0:
                        print(f"  Processed {i:,} fragments...")
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        chrom, start, end, barcode, count = parts[:5]
                        fragments.append({
                            'chrom': chrom,
                            'start': int(start),
                            'end': int(end),
                            'barcode': barcode,
                            'count': int(count) if len(parts) > 4 else 1
                        })
                        
                        # Count fragments per barcode
                        if barcode not in barcode_counts:
                            barcode_counts[barcode] = 0
                        barcode_counts[barcode] += int(count) if len(parts) > 4 else 1
                    
                    # Limit initial loading for memory efficiency
                    if i > 2000000:  # Load first 2M fragments for analysis
                        break
            
            fragments_df = pd.DataFrame(fragments)
            barcode_stats = pd.DataFrame([
                {'barcode': bc, 'fragment_count': count} 
                for bc, count in barcode_counts.items()
            ]).sort_values('fragment_count', ascending=False)
            
            print(f"  Loaded {len(fragments):,} fragments from {len(barcode_counts):,} unique barcodes")
            print(f"  Sample barcodes: {barcode_stats['barcode'].head().tolist()}")
            print(f"  Fragment counts per barcode: {barcode_stats['fragment_count'].describe()}")
            
            return fragments_df, barcode_stats
            
        except Exception as e:
            print(f"Error loading ATAC data for {sample_name}: {e}")
            return None, None
    
    def analyze_barcode_overlap(self, rna_sample, atac_sample, rna_adata, atac_stats):
        """Analyze overlap between RNA and ATAC barcodes"""
        print(f"\nAnalyzing barcode overlap between RNA and ATAC for sample comparison...")
        
        # Get barcode sets
        rna_barcodes = set(rna_adata.obs['barcode_clean'])
        atac_barcodes = set(atac_stats['barcode'])
        
        # Calculate overlaps
        overlap = rna_barcodes.intersection(atac_barcodes)
        rna_only = rna_barcodes - atac_barcodes
        atac_only = atac_barcodes - rna_barcodes
        
        # Create summary
        summary = {
            'rna_sample': rna_sample,
            'atac_sample': atac_sample,
            'rna_total': len(rna_barcodes),
            'atac_total': len(atac_barcodes),
            'overlap_count': len(overlap),
            'rna_only_count': len(rna_only),
            'atac_only_count': len(atac_only),
            'overlap_fraction_rna': len(overlap) / len(rna_barcodes) if rna_barcodes else 0,
            'overlap_fraction_atac': len(overlap) / len(atac_barcodes) if atac_barcodes else 0
        }
        
        print(f"Overlap Analysis Results:")
        print(f"  RNA cells: {summary['rna_total']:,}")
        print(f"  ATAC cells: {summary['atac_total']:,}")
        print(f"  Overlapping cells: {summary['overlap_count']:,}")
        print(f"  RNA overlap fraction: {summary['overlap_fraction_rna']:.3f}")
        print(f"  ATAC overlap fraction: {summary['overlap_fraction_atac']:.3f}")
        
        return summary, overlap, rna_only, atac_only
    
    def create_coupled_dataset(self, rna_adata, atac_stats, overlap_barcodes, sample_name):
        """Create coupled dataset with both RNA and ATAC information"""
        print(f"\nCreating coupled dataset for {sample_name}...")
        
        # Filter RNA data to overlapping barcodes
        rna_overlap_mask = rna_adata.obs['barcode_clean'].isin(overlap_barcodes)
        rna_coupled = rna_adata[rna_overlap_mask, :].copy()
        
        # Filter ATAC data to overlapping barcodes
        atac_overlap = atac_stats[atac_stats['barcode'].isin(overlap_barcodes)].copy()
        
        # Merge ATAC information into RNA obs
        atac_dict = atac_overlap.set_index('barcode')['fragment_count'].to_dict()
        rna_coupled.obs['atac_fragments'] = rna_coupled.obs['barcode_clean'].map(atac_dict)
        rna_coupled.obs['has_atac'] = True
        rna_coupled.obs['coupled_sample'] = sample_name
        
        print(f"  Coupled dataset: {rna_coupled.n_obs} cells x {rna_coupled.n_vars} genes")
        print(f"  ATAC fragment distribution: {rna_coupled.obs['atac_fragments'].describe()}")
        
        return rna_coupled
    
    def run_coupling_analysis(self):
        """Run complete coupling analysis pipeline"""
        print("="*80)
        print("MULTIMODAL COUPLING ANALYSIS")
        print("="*80)
        
        # Find data paths
        rna_paths = self.find_rna_data_paths()
        atac_paths = self.find_atac_data_paths()
        
        print(f"\nFound RNA datasets: {list(rna_paths.keys())}")
        print(f"Found ATAC datasets: {list(atac_paths.keys())}")
        
        # Load all datasets
        print("\n" + "="*50)
        print("LOADING DATA")
        print("="*50)
        
        for sample, path in rna_paths.items():
            adata = self.load_rna_data(sample, path)
            if adata is not None:
                self.rna_data[sample] = adata
        
        for sample, path in atac_paths.items():
            fragments_df, barcode_stats = self.load_atac_fragments(sample, path)
            if fragments_df is not None:
                self.atac_data[sample] = {
                    'fragments': fragments_df,
                    'barcode_stats': barcode_stats
                }
        
        # Analyze overlaps
        print("\n" + "="*50)
        print("BARCODE OVERLAP ANALYSIS")
        print("="*50)
        
        overlap_results = []
        
        for rna_sample in self.rna_data.keys():
            for atac_sample in self.atac_data.keys():
                # Try to match samples (accounting for naming differences)
                if self._samples_match(rna_sample, atac_sample):
                    print(f"\nMatching {rna_sample} (RNA) with {atac_sample} (ATAC)")
                    
                    summary, overlap, rna_only, atac_only = self.analyze_barcode_overlap(
                        rna_sample, atac_sample, 
                        self.rna_data[rna_sample], 
                        self.atac_data[atac_sample]['barcode_stats']
                    )
                    
                    overlap_results.append(summary)
                    
                    # Create coupled dataset if significant overlap
                    if len(overlap) > 100:  # Minimum threshold
                        coupled_name = f"{rna_sample}_coupled"
                        coupled_data = self.create_coupled_dataset(
                            self.rna_data[rna_sample],
                            self.atac_data[atac_sample]['barcode_stats'],
                            overlap,
                            coupled_name
                        )
                        self.coupled_data[coupled_name] = coupled_data
        
        # Store results
        self.results['overlap_analysis'] = pd.DataFrame(overlap_results)
        
        return self.results
    
    def _samples_match(self, rna_sample, atac_sample):
        """Check if RNA and ATAC samples correspond to the same biological sample"""
        # Simple matching based on sample name patterns
        rna_clean = rna_sample.replace('_', '-').lower()
        atac_clean = atac_sample.replace('_', '-').lower()
        
        # Extract key identifiers
        rna_parts = set(rna_clean.split('-'))
        atac_parts = set(atac_clean.split('-'))
        
        # Check for significant overlap in name parts
        common_parts = rna_parts.intersection(atac_parts)
        return len(common_parts) >= 3  # At least 3 matching parts (e.g., r26, nestin, ctrl/mut)
    
    def generate_qc_plots(self, output_dir=None):
        """Generate quality control plots for coupling analysis"""
        if output_dir is None:
            output_dir = self.base_dir / "coupling_analysis_output"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        print(f"\nGenerating QC plots in: {output_dir}")
        
        # Plot 1: Barcode overlap summary
        if 'overlap_analysis' in self.results:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            
            df = self.results['overlap_analysis']
            
            # Overlap counts
            ax1.bar(range(len(df)), df['overlap_count'])
            ax1.set_xlabel('Sample Pairs')
            ax1.set_ylabel('Overlapping Cells')
            ax1.set_title('Number of Overlapping Cells')
            ax1.set_xticks(range(len(df)))
            ax1.set_xticklabels([f"{r}\nvs\n{a}" for r, a in zip(df['rna_sample'], df['atac_sample'])], 
                               rotation=0, ha='center')
            
            # Overlap fractions
            x = range(len(df))
            ax2.bar([i-0.2 for i in x], df['overlap_fraction_rna'], width=0.4, label='RNA fraction', alpha=0.7)
            ax2.bar([i+0.2 for i in x], df['overlap_fraction_atac'], width=0.4, label='ATAC fraction', alpha=0.7)
            ax2.set_xlabel('Sample Pairs')
            ax2.set_ylabel('Overlap Fraction')
            ax2.set_title('Fraction of Overlapping Cells')
            ax2.set_xticks(x)
            ax2.set_xticklabels([f"{r}\nvs\n{a}" for r, a in zip(df['rna_sample'], df['atac_sample'])], 
                               rotation=0, ha='center')
            ax2.legend()
            
            plt.tight_layout()
            plt.savefig(output_dir / "barcode_overlap_summary.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        # Plot 2: Individual sample QC
        for name, adata in self.coupled_data.items():
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            
            # RNA QC metrics
            adata.var['mt'] = adata.var_names.str.startswith('Mt-')
            sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
            
            # Total counts
            axes[0,0].hist(adata.obs['total_counts'], bins=50, alpha=0.7)
            axes[0,0].set_xlabel('Total RNA counts')
            axes[0,0].set_ylabel('Cells')
            axes[0,0].set_title('RNA: Total Counts Distribution')
            
            # Number of genes
            axes[0,1].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7)
            axes[0,1].set_xlabel('Number of genes')
            axes[0,1].set_ylabel('Cells')
            axes[0,1].set_title('RNA: Gene Count Distribution')
            
            # ATAC fragments
            axes[1,0].hist(adata.obs['atac_fragments'], bins=50, alpha=0.7)
            axes[1,0].set_xlabel('ATAC fragments')
            axes[1,0].set_ylabel('Cells')
            axes[1,0].set_title('ATAC: Fragment Count Distribution')
            
            # Correlation
            axes[1,1].scatter(adata.obs['total_counts'], adata.obs['atac_fragments'], alpha=0.5)
            axes[1,1].set_xlabel('RNA total counts')
            axes[1,1].set_ylabel('ATAC fragments')
            axes[1,1].set_title('RNA vs ATAC Correlation')
            
            plt.suptitle(f'Quality Control: {name}')
            plt.tight_layout()
            plt.savefig(output_dir / f"qc_{name}.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        print(f"QC plots saved to: {output_dir}")
    
    def save_coupled_datasets(self, output_dir=None):
        """Save coupled datasets for downstream analysis"""
        if output_dir is None:
            output_dir = self.base_dir / "coupling_analysis_output"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        print(f"\nSaving coupled datasets to: {output_dir}")
        
        for name, adata in self.coupled_data.items():
            # Save as h5ad
            output_path = output_dir / f"{name}.h5ad"
            adata.write(output_path)
            print(f"  Saved {name}: {output_path}")
            
            # Also save barcode lists
            barcode_path = output_dir / f"{name}_barcodes.txt"
            with open(barcode_path, 'w') as f:
                for barcode in adata.obs['barcode_clean']:
                    f.write(f"{barcode}\n")
            print(f"  Saved barcodes: {barcode_path}")
        
        # Save overlap analysis results
        if 'overlap_analysis' in self.results:
            results_path = output_dir / "overlap_analysis_results.csv"
            self.results['overlap_analysis'].to_csv(results_path, index=False)
            print(f"  Saved overlap analysis: {results_path}")

def main():
    """Main execution function"""
    print("Starting Multimodal Coupling Analysis...")
    
    # Initialize coupler
    coupler = MultimodalCoupler()
    
    # Run analysis
    results = coupler.run_coupling_analysis()
    
    # Generate outputs
    coupler.generate_qc_plots()
    coupler.save_coupled_datasets()
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    
    if 'overlap_analysis' in results:
        print("\nSUMMARY:")
        for _, row in results['overlap_analysis'].iterrows():
            print(f"Sample pair: {row['rna_sample']} vs {row['atac_sample']}")
            print(f"  Overlapping cells: {row['overlap_count']:,}")
            print(f"  RNA overlap: {row['overlap_fraction_rna']:.1%}")
            print(f"  ATAC overlap: {row['overlap_fraction_atac']:.1%}")
    
    print(f"\nCoupled datasets created: {len(coupler.coupled_data)}")
    for name in coupler.coupled_data.keys():
        print(f"  - {name}")
    
    print("\nNext steps:")
    print("1. Review QC plots in coupling_analysis_output/")
    print("2. Load coupled datasets (.h5ad files) for downstream analysis")
    print("3. Run integration analysis using Seurat/Scanpy")

if __name__ == "__main__":
    main()