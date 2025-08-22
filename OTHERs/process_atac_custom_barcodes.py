#!/usr/bin/env python3
"""
Custom ATAC-seq barcode processing script
Handles ARC-v1 multiome data where barcodes are in R2 instead of R1
"""

import argparse
import gzip
import os
from collections import defaultdict
import sys

def read_whitelist(whitelist_file):
    """Read barcode whitelist from file"""
    barcodes = set()
    with open(whitelist_file, 'r') as f:
        for line in f:
            barcode = line.strip()
            if barcode:
                barcodes.add(barcode)
    return barcodes

def analyze_barcodes(sample, data_dir, whitelist_file, output_dir):
    """Analyze barcode positions in ATAC R2 reads"""
    print(f"Analyzing barcode positions for {sample}")
    
    # Read whitelist
    valid_barcodes = read_whitelist(whitelist_file)
    print(f"Loaded {len(valid_barcodes)} valid barcodes from whitelist")
    
    # Files to analyze
    r2_file = os.path.join(data_dir, f"{sample}_R2_001.fastq.gz")
    
    if not os.path.exists(r2_file):
        print(f"ERROR: R2 file not found: {r2_file}")
        return False
    
    # Analyze first 100000 reads
    barcode_positions = defaultdict(int)
    total_reads = 0
    
    print("Analyzing R2 reads for barcode positions...")
    
    with gzip.open(r2_file, 'rt') as f:
        while total_reads < 100000:
            try:
                header = f.readline()
                if not header:
                    break
                sequence = f.readline().strip()
                plus = f.readline()
                quality = f.readline()
                
                total_reads += 1
                
                # Check different barcode positions in the read
                for start_pos in range(0, min(len(sequence) - 10, 16)):
                    for length in [10, 12, 16]:
                        if start_pos + length <= len(sequence):
                            potential_barcode = sequence[start_pos:start_pos + length]
                            if potential_barcode in valid_barcodes:
                                position_key = f"pos_{start_pos}_len_{length}"
                                barcode_positions[position_key] += 1
                
            except Exception as e:
                print(f"Error reading read {total_reads}: {e}")
                break
    
    # Report results
    print(f"\nAnalyzed {total_reads} reads")
    print("Barcode position analysis:")
    
    if not barcode_positions:
        print("No valid barcodes found in any position")
        return False
    
    # Find best position
    best_position = max(barcode_positions.keys(), key=lambda x: barcode_positions[x])
    best_count = barcode_positions[best_position]
    
    for position, count in sorted(barcode_positions.items(), key=lambda x: x[1], reverse=True):
        percentage = (count / total_reads) * 100
        print(f"  {position}: {count} matches ({percentage:.2f}%)")
    
    print(f"\nBest position: {best_position} with {best_count} matches")
    
    # Save analysis results
    analysis_file = os.path.join(output_dir, f"{sample}_barcode_analysis.txt")
    with open(analysis_file, 'w') as f:
        f.write(f"Sample: {sample}\n")
        f.write(f"Total reads analyzed: {total_reads}\n")
        f.write(f"Best barcode position: {best_position}\n")
        f.write(f"Best position matches: {best_count}\n")
        f.write("\nAll positions:\n")
        for position, count in sorted(barcode_positions.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_reads) * 100
            f.write(f"{position}: {count} matches ({percentage:.2f}%)\n")
    
    print(f"Analysis saved to: {analysis_file}")
    return True

def filter_reads(sample, data_dir, whitelist_file, output_dir):
    """Filter ATAC reads using valid RNA barcodes"""
    print(f"Filtering reads for {sample}")
    
    # Read whitelist
    valid_barcodes = read_whitelist(whitelist_file)
    print(f"Using {len(valid_barcodes)} valid barcodes for filtering")
    
    # Input files
    files = {
        'I1': os.path.join(data_dir, f"{sample}_I1_001.fastq.gz"),
        'R1': os.path.join(data_dir, f"{sample}_R1_001.fastq.gz"),
        'R2': os.path.join(data_dir, f"{sample}_R2_001.fastq.gz"),
        'R3': os.path.join(data_dir, f"{sample}_R3_001.fastq.gz")
    }
    
    # Check input files exist
    for read_type, filepath in files.items():
        if not os.path.exists(filepath):
            print(f"ERROR: Input file not found: {filepath}")
            return False
    
    # Output files
    output_files = {
        'I1': os.path.join(output_dir, f"{sample}_filtered_I1_001.fastq.gz"),
        'R1': os.path.join(output_dir, f"{sample}_filtered_R1_001.fastq.gz"),
        'R2': os.path.join(output_dir, f"{sample}_filtered_R2_001.fastq.gz"),
        'R3': os.path.join(output_dir, f"{sample}_filtered_R3_001.fastq.gz")
    }
    
    print("Opening input and output files...")
    
    # Open all files
    input_handles = {}
    output_handles = {}
    
    try:
        for read_type in ['I1', 'R1', 'R2', 'R3']:
            input_handles[read_type] = gzip.open(files[read_type], 'rt')
            output_handles[read_type] = gzip.open(output_files[read_type], 'wt')
        
        total_reads = 0
        filtered_reads = 0
        
        print("Processing reads...")
        
        while True:
            # Read one record from each file
            records = {}
            all_good = True
            
            for read_type in ['I1', 'R1', 'R2', 'R3']:
                try:
                    header = input_handles[read_type].readline()
                    if not header:
                        all_good = False
                        break
                    sequence = input_handles[read_type].readline()
                    plus = input_handles[read_type].readline()
                    quality = input_handles[read_type].readline()
                    
                    records[read_type] = (header, sequence, plus, quality)
                except:
                    all_good = False
                    break
            
            if not all_good:
                break
            
            total_reads += 1
            
            # Extract barcode from R2 (assuming position 0, length 16)
            r2_sequence = records['R2'][1].strip()
            if len(r2_sequence) >= 16:
                barcode = r2_sequence[:16]
                
                # Check if barcode is valid (try different lengths)
                barcode_found = False
                for length in [16, 12, 10]:
                    test_barcode = r2_sequence[:length]
                    if test_barcode in valid_barcodes:
                        barcode_found = True
                        break
                
                if barcode_found:
                    # Write filtered records
                    for read_type in ['I1', 'R1', 'R2', 'R3']:
                        for line in records[read_type]:
                            output_handles[read_type].write(line)
                    filtered_reads += 1
            
            if total_reads % 1000000 == 0:
                print(f"Processed {total_reads} reads, kept {filtered_reads}")
    
    finally:
        # Close all files
        for handle in input_handles.values():
            handle.close()
        for handle in output_handles.values():
            handle.close()
    
    print(f"\nFiltering complete:")
    print(f"Total reads processed: {total_reads}")
    print(f"Reads with valid barcodes: {filtered_reads}")
    print(f"Filtering rate: {(filtered_reads/total_reads)*100:.2f}%")
    
    # Save filtering stats
    stats_file = os.path.join(output_dir, f"{sample}_filtering_stats.txt")
    with open(stats_file, 'w') as f:
        f.write(f"Sample: {sample}\n")
        f.write(f"Total reads: {total_reads}\n")
        f.write(f"Filtered reads: {filtered_reads}\n")
        f.write(f"Filtering rate: {(filtered_reads/total_reads)*100:.2f}%\n")
    
    print(f"Stats saved to: {stats_file}")
    return True

def main():
    parser = argparse.ArgumentParser(description='Process ATAC data with custom barcode handling')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--data-dir', required=True, help='Directory containing FASTQ files')
    parser.add_argument('--whitelist', required=True, help='Barcode whitelist file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--analyze-only', action='store_true', help='Only analyze, do not filter')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.analyze_only:
        success = analyze_barcodes(args.sample, args.data_dir, args.whitelist, args.output_dir)
    else:
        success = filter_reads(args.sample, args.data_dir, args.whitelist, args.output_dir)
    
    if not success:
        sys.exit(1)
    
    print("Processing completed successfully")

if __name__ == '__main__':
    main()