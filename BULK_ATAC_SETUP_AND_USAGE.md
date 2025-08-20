# Bulk ATAC-seq Processing - Setup and Usage Guide

**Date**: August 20, 2025  
**Purpose**: Process ATAC-seq data as bulk samples, bypassing single-cell barcode issues  
**Script**: [`run_bulk_atac_processing.sh`](ATAC_data/run_bulk_atac_processing.sh)

---

## 🎯 Overview

This approach treats your ATAC-seq data as **bulk samples** rather than single-cell, which completely bypasses the barcode detection problems that prevented single-cell analysis.

**What it does**:
- Uses genomic reads from R1 (50bp) and R3 (49bp) 
- Ignores R2 (barcode reads) entirely
- Performs standard bulk ATAC-seq analysis
- Compares Ctrl vs Mut at population level

---

## 📋 Prerequisites Check

Before running, verify these tools are available:

### Required Tools:
```bash
# Check if tools are available
which bowtie2      # For alignment
which samtools     # For BAM processing  
which macs2        # For peak calling

# Check Java (needed for Trimmomatic/Picard)
java -version
```

### Optional Tools (script will work without them):
- **Trimmomatic**: For quality trimming (script will skip if not found)
- **Picard**: For duplicate removal (script will skip if not found)

---

## 🚀 Quick Start

### 1. **Submit the Job**
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data
sbatch run_bulk_atac_processing.sh
```

### 2. **Monitor Progress**
```bash
# Check job status
squeue -u kubacki.michal

# Monitor logs
tail -f logs/bulk_atac_seq_*.out
tail -f logs/bulk_atac_seq_*.err
```

---

## 📁 Expected Output Structure

```
bulk_atac_output/
├── R26-Nestin-Ctrl-adult/
│   ├── trimmed/                    # Quality-trimmed reads
│   ├── alignment/
│   │   ├── R26-Nestin-Ctrl-adult_sorted.bam
│   │   ├── R26-Nestin-Ctrl-adult_dedup.bam    # Final BAM file
│   │   └── *.bam.bai                           # Index files
│   ├── peaks/
│   │   ├── R26-Nestin-Ctrl-adult_peaks.narrowPeak  # Peak calls
│   │   ├── R26-Nestin-Ctrl-adult_peaks.broadPeak
│   │   ├── R26-Nestin-Ctrl-adult_summits.bed
│   │   └── R26-Nestin-Ctrl-adult_model.r
│   ├── stats/
│   │   ├── R26-Nestin-Ctrl-adult_flagstat.txt      # Alignment stats
│   │   └── R26-Nestin-Ctrl-adult_idxstats.txt
│   ├── tracks/                     # Browser tracks (optional)
│   └── R26-Nestin-Ctrl-adult_processing_report.txt
└── R26-Nestin-Mut-adult/
    └── [same structure]
```

---

## ✅ Success Indicators

### 1. **Check Processing Reports**
```bash
cat bulk_atac_output/R26-Nestin-Ctrl-adult/R26-Nestin-Ctrl-adult_processing_report.txt
cat bulk_atac_output/R26-Nestin-Mut-adult/R26-Nestin-Mut-adult_processing_report.txt
```

### 2. **Expected Results**
- **Alignment rate**: >70% (check flagstat files)
- **Peak count**: 20,000-100,000 peaks per sample
- **BAM files**: >1GB each
- **Processing reports**: Complete statistics

### 3. **Quick Quality Check**
```bash
# Check number of peaks called
wc -l bulk_atac_output/*/peaks/*_peaks.narrowPeak

# Check alignment statistics
grep "mapped (" bulk_atac_output/*/stats/*_flagstat.txt
```

---

## 🔧 Customization Options

### Modify Tool Paths
Edit these variables in the script if tools are in different locations:
```bash
TRIMMOMATIC="/path/to/trimmomatic.jar"
BOWTIE2="/path/to/bowtie2"
PICARD="/path/to/picard.jar"
```

### Adjust MACS2 Parameters
For different peak calling sensitivity:
```bash
# More stringent (fewer peaks)
--broad-cutoff 0.05

# More permissive (more peaks)  
--broad-cutoff 0.2
```

---

## 📊 Downstream Analysis

### 1. **Differential Peak Analysis**
Use tools like:
- **DiffBind** (R package) - Compare peaks between Ctrl vs Mut
- **DESeq2** - Statistical analysis of count differences
- **ChIPseeker** - Peak annotation and visualization

### 2. **Integration with RNA Data**
- Load both ATAC peaks and RNA results into R/Python
- Correlate accessibility changes with expression changes
- Identify regulatory relationships

### 3. **Visualization**
- Load BAM files into **IGV** browser
- Create track hub for **UCSC Genome Browser**  
- Generate heatmaps of peak signals

---

## 🚨 Troubleshooting

### Common Issues:

#### **Job fails immediately**
```bash
# Check tool availability
module load samtools bowtie2 # if using modules
which macs2
```

#### **Low alignment rates (<50%)**
- Check if Bowtie2 index path is correct
- Verify reference genome version matches your samples

#### **No peaks called**
- Check MACS2 output in logs
- Try less stringent parameters (`--broad-cutoff 0.2`)

#### **Out of memory errors**
- Increase memory allocation: `#SBATCH --mem=128G`
- Process samples individually instead of array job

---

## 💡 Why This Approach Works

**Bypasses Single-Cell Issues**:
- ✅ No barcode detection required
- ✅ Uses actual genomic DNA reads (R1 + R3)
- ✅ Standard, well-established workflow
- ✅ Compatible with downstream analysis tools

**Biological Insights Still Possible**:
- ✅ Differential accessibility between Ctrl vs Mut
- ✅ Peak annotation and gene association  
- ✅ Integration with RNA expression data
- ✅ Pathway and GO term analysis

---

## 🎯 Expected Timeline

- **Per sample processing**: 2-6 hours
- **Both samples (parallel)**: 4-8 hours total
- **Peak calling**: 30 minutes - 2 hours
- **Total runtime**: 6-10 hours

---

## 🏁 Next Steps After Processing

1. **Verify results** using success indicators above
2. **Compare peak counts** between Ctrl and Mut samples  
3. **Load data into R/Python** for differential analysis
4. **Integrate with RNA results** for comprehensive analysis
5. **Generate publication-quality figures**

---

**📧 Status**: Ready to run - bulk processing approach bypasses all single-cell issues  
**🎯 Expected Outcome**: High-quality bulk ATAC-seq results suitable for differential analysis