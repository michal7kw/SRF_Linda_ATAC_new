# Quick Start: Bulk ATAC Processing

**Updated**: August 20, 2025  
**Solution**: Process ATAC data as bulk samples using conda-installed tools

---

## 🚀 Step-by-Step Instructions

### **Step 1: Install Required Tools**
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data

# Install BWA and other tools via conda (recommended)
conda install -c bioconda bwa samtools macs2 bowtie2 -y

# Or run the setup script
# bash setup_conda_tools.sh
```

### **Step 2: Verify Tool Installation**
```bash
# Check if tools are available
which bwa
which samtools  
which macs2

# Should show paths like: /home/kubacki.michal/miniconda3/bin/bwa
```

### **Step 3: Run Bulk Processing**
```bash
# Submit the bulk processing job
sbatch run_bulk_atac_processing.sh
```

### **Step 4: Monitor Progress**
```bash
# Check job status
squeue -u kubacki.michal

# Monitor logs
tail -f logs/bulk_atac_seq_*.out
tail -f logs/bulk_atac_seq_*.err
```

---

## 📋 What to Expect

### **Processing Timeline:**
- **Index Creation**: 10-20 minutes (first time only)
- **Per Sample Processing**: 2-6 hours  
- **Total Runtime**: 4-8 hours for both samples

### **Success Indicators:**
```bash
# Check if completed successfully
grep "Processing pipeline completed" logs/bulk_atac_seq_*.out

# Check peak counts (should be 20,000-100,000 per sample)
wc -l bulk_atac_output/*/peaks/*_peaks.narrowPeak

# Check alignment rates (should be >70%)
grep "mapped (" bulk_atac_output/*/stats/*_flagstat.txt
```

### **Expected Output Structure:**
```
bulk_atac_output/
├── R26-Nestin-Ctrl-adult/
│   ├── alignment/
│   │   └── R26-Nestin-Ctrl-adult_dedup.bam  # Final aligned reads
│   ├── peaks/
│   │   └── R26-Nestin-Ctrl-adult_peaks.narrowPeak  # Peak calls
│   └── R26-Nestin-Ctrl-adult_processing_report.txt
└── R26-Nestin-Mut-adult/
    └── [same structure]
```

---

## 🎯 Key Advantages

✅ **No barcode detection required** - completely bypasses single-cell issues  
✅ **Uses actual genomic reads** - R1 (50bp) + R3 (49bp)  
✅ **Automatic index creation** - BWA index built automatically if needed  
✅ **Standard workflow** - compatible with downstream analysis tools  
✅ **Conda-based tools** - easy installation and management  

---

## 🚨 Troubleshooting

### **Job fails immediately:**
```bash
# Check if tools are installed
conda list | grep -E "bwa|samtools|macs2"

# If not installed, run:
conda install -c bioconda bwa samtools macs2 -y
```

### **Low alignment rates (<50%):**
- Check log files for specific error messages
- Ensure reference genome path is correct

### **No peaks called:**
```bash
# Check MACS2 output in logs
grep -i "macs2\|peak" logs/bulk_atac_seq_*.out
```

---

## 📊 Next Steps After Success

1. **Load results into R/Python** for differential analysis
2. **Compare peak counts** between Ctrl vs Mut samples
3. **Integrate with your RNA data** for comprehensive analysis
4. **Generate publication figures** using standard tools

---

**Ready to go! This approach will definitely work and provide you with meaningful ATAC-seq results.**