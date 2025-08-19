#!/usr/bin/env python3
"""
Custom ATAC processing for separately sequenced ARC-v1 data
Extracts barcodes from R2 and processes using RNA whitelist
"""

import os
import sys
import gzip
import argparse
from pathlib import Path
from collections import Counter
import subprocess

def extract_barcode_from_r2(r2_file, start=0, length=16, max_reads=1000000):
    """
    Extract potential cell barcodes from R2 reads
    In ARC-v1, barcodes might be in the first 16bp of R2
    """
    barcodes = []
    with gzip.open(r2_file, 'rt') as f:
        line_num = 0
        for line in f:
            if line_num % 4 == 1:  # Sequence line
                barcode = line[start:start+length].strip()
                if len(barcode) == length:  # Valid length
                    barcodes.append(barcode)
            line_num += 1
            if len(barcodes) >= max_reads:
                break
    return barcodes

def load_whitelist(whitelist_file):
    """Load barcode whitelist from RNA processing"""
    valid_barcodes = set()
    with open(whitelist_file, 'r') as f:
        for line in f:
            barcode = line.strip().split('-')[0]  # Remove -1 suffix if present
            if barcode:
                valid_barcodes.add(barcode)
    return valid_barcodes

def analyze_barcode_positions(r2_file, whitelist, max_reads=100000):
    """
    Try different positions in R2 to find where barcodes are located
    """
    print("Analyzing different barcode positions in R2...")
    
    # Read sample of R2
    sequences = []
    with gzip.open(r2_file, 'rt') as f:
        line_num = 0
        for line in f:
            if line_num % 4 == 1:  # Sequence line
                sequences.append(line.strip())
            line_num += 1
            if len(sequences) >= max_reads:
                break
    
    # Try different positions and lengths
    results = []
    for start in [0, 8]:  # Try start at 0 or 8
        for length in [10, 12, 16]:  # Common barcode lengths
            barcodes = [seq[start:start+length] for seq in sequences if len(seq) >= start+length]
            if barcodes:
                # Count matches with whitelist
                matches = sum(1 for bc in barcodes if bc in whitelist)
                match_rate = matches / len(barcodes) * 100
                results.append({
                    'start': start,
                    'length': length,
                    'matches': matches,
                    'total': len(barcodes),
                    'match_rate': match_rate
                })
                print(f"  Position {start}:{start+length} - {matches}/{len(barcodes)} ({match_rate:.1f}%) match RNA barcodes")
    
    # Find best position
    if results:
        best = max(results, key=lambda x: x['match_rate'])
        print(f"\nBest barcode position: {best['start']}:{best['start']+best['length']} with {best['match_rate']:.1f}% match rate")
        return best['start'], best['length']
    return 0, 16

def create_filtered_fastq(r1_file, r2_file, r3_file, i1_file, whitelist, output_dir, 
                         barcode_start=0, barcode_length=16):
    """
    Create new FASTQ files filtered by valid barcodes from RNA
    """
    sample_name = Path(r1_file).stem.replace('_R1_001.fastq', '')
    
    # Output files
    out_r1 = output_dir / f"{sample_name}_filtered_R1_001.fastq.gz"
    out_r2 = output_dir / f"{sample_name}_filtered_R2_001.fastq.gz"
    out_r3 = output_dir / f"{sample_name}_filtered_R3_001.fastq.gz"
    out_i1 = output_dir / f"{sample_name}_filtered_I1_001.fastq.gz"
    
    total_reads = 0
    valid_reads = 0
    
    print(f"Filtering reads with valid barcodes...")
    
    with gzip.open(r1_file, 'rt') as f1, \
         gzip.open(r2_file, 'rt') as f2, \
         gzip.open(r3_file, 'rt') as f3, \
         gzip.open(i1_file, 'rt') as fi, \
         gzip.open(out_r1, 'wt') as o1, \
         gzip.open(out_r2, 'wt') as o2, \
         gzip.open(out_r3, 'wt') as o3, \
         gzip.open(out_i1, 'wt') as oi:
        
        while True:
            # Read 4 lines from each file (one read)
            r1_lines = [f1.readline() for _ in range(4)]
            r2_lines = [f2.readline() for _ in range(4)]
            r3_lines = [f3.readline() for _ in range(4)]
            i1_lines = [fi.readline() for _ in range(4)]
            
            if not r1_lines[0]:  # End of file
                break
            
            total_reads += 1
            
            # Extract barcode from R2
            barcode = r2_lines[1][barcode_start:barcode_start+barcode_length].strip()
            
            # Check if barcode is in whitelist
            if barcode in whitelist:
                valid_reads += 1
                # Write all reads
                for line in r1_lines:
                    o1.write(line)
                for line in r2_lines:
                    o2.write(line)
                for line in r3_lines:
                    o3.write(line)
                for line in i1_lines:
                    oi.write(line)
            
            if total_reads % 100000 == 0:
                print(f"  Processed {total_reads} reads, kept {valid_reads} ({valid_reads/total_reads*100:.1f}%)")
    
    print(f"Filtering complete: {valid_reads}/{total_reads} ({valid_reads/total_reads*100:.1f}%) reads kept")
    return valid_reads, total_reads

def main():
    parser = argparse.ArgumentParser(description='Process ATAC data with RNA barcodes')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--data-dir', required=True, help='Directory with ATAC FASTQ files')
    parser.add_argument('--whitelist', required=True, help='RNA barcode whitelist file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--analyze-only', action='store_true', help='Only analyze barcode positions')
    
    args = parser.parse_args()
    
    # Setup paths
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find input files
    r1_file = data_dir / f"{args.sample}_R1_001.fastq.gz"
    r2_file = data_dir / f"{args.sample}_R2_001.fastq.gz"
    r3_file = data_dir / f"{args.sample}_R3_001.fastq.gz"
    i1_file = data_dir / f"{args.sample}_I1_001.fastq.gz"
    
    # Check files exist
    for f in [r1_file, r2_file, r3_file, i1_file]:
        if not f.exists():
            print(f"ERROR: File not found: {f}")
            sys.exit(1)
    
    print(f"Processing sample: {args.sample}")
    print(f"Using RNA whitelist: {args.whitelist}")
    
    # Load whitelist
    whitelist = load_whitelist(args.whitelist)
    print(f"Loaded {len(whitelist)} valid barcodes from RNA")
    
    # Analyze barcode positions
    barcode_start, barcode_length = analyze_barcode_positions(r2_file, whitelist)
    
    if not args.analyze_only:
        # Filter reads
        valid, total = create_filtered_fastq(
            r1_file, r2_file, r3_file, i1_file,
            whitelist, output_dir,
            barcode_start, barcode_length
        )
        
        print(f"\nResults saved to: {output_dir}")
        print(f"You can now process these filtered files with standard ATAC-seq pipelines")

if __name__ == "__main__":
    main()