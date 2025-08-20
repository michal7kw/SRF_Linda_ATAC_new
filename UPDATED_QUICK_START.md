# Updated Quick Start: Bulk ATAC Processing
**Fixed**: Conda dependency conflicts resolved  
**Date**: August 20, 2025

---

## 🛠️ Tool Installation (Fixed for Python 3.12)

### **Method 1: Step-by-step Installation (Recommended)**
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data

# Install BWA and samtools via conda
conda install -c bioconda bwa samtools -y

# Install MACS2 via pip (avoids conda conflicts)
pip install MACS2

# Verify installations
which bwa samtools macs2
```

### **Method 2: Use Alternative Installation Script**
```bash
bash install_tools_alternative.sh
```

### **Method 3: Manual Installation (if above fails)**
```bash
# Just install BWA and samtools (minimum requirements)
conda install -c bioconda bwa samtools -y

# MACS2 is optional - we can do basic peak calling without it
```

---

## 🚀 Run Bulk Processing

```bash
# After tools are installed, submit the job
sbatch run_bulk_atac_processing.sh
```

---

## 📋 Monitor and Check Results

```bash
# Monitor job
squeue -u kubacki.michal
tail -f logs/bulk_atac_seq_*.out

# Check success after completion
grep "Processing pipeline completed" logs/bulk_atac_seq_*.out
ls -la bulk_atac_output/*/alignment/*.bam
ls -la bulk_atac_output/*/peaks/
```

---

## 🎯 What Makes This Work Now

✅ **BWA + samtools**: Core alignment tools installed via conda  
✅ **MACS2 via pip**: Avoids Python 3.12 conflicts  
✅ **Automatic indexing**: Script creates BWA index automatically  
✅ **No barcode issues**: Uses genomic reads R1+R3 only  
✅ **Robust error handling**: Works even if some tools are missing  

---

## 💡 Alternative if Tools Still Fail

If you're still having installation issues, the script can work with just BWA and samtools:

```bash
# Minimal installation
conda install -c bioconda bwa samtools -y

# The script will automatically:
# - Use BWA for alignment
# - Skip MACS2 peak calling if not available
# - Still generate useful BAM files for manual analysis
```

You can then do peak calling manually later:
```bash
# After alignment completes, you can install MACS2 separately
pip install MACS2

# Then run peak calling manually on the BAM files
macs2 callpeak -t bulk_atac_output/R26-Nestin-Ctrl-adult/alignment/*.bam \
               -f BAMPE -g mm -n Ctrl --outdir peaks_manual/
```

---

**The key point is that BWA + samtools are sufficient to generate aligned BAM files, and peak calling can be done separately if needed.**