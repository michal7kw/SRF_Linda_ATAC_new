# Final Setup Guide: Bulk ATAC Processing with Conda Environment

**Date**: August 20, 2025  
**Solution**: Dedicated conda environment with Python 3.10 to avoid dependency conflicts

---

## 🎯 Complete Setup and Usage

### **Step 1: Create Bioinformatics Environment**
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data

# Run the setup script (creates environment and installs all tools)
bash setup_bioinformatics_env.sh
```

**What this does:**
- Creates new conda environment `atac_analysis` with Python 3.10
- Installs BWA, samtools, MACS2, bowtie2 in the new environment
- Avoids conflicts with your base Python 3.12 environment
- Verifies all installations

### **Step 2: Run Bulk Processing**
```bash
# The script automatically activates the conda environment
sbatch run_bulk_atac_processing.sh
```

### **Step 3: Monitor Progress**
```bash
# Check job status
squeue -u kubacki.michal

# Monitor logs
tail -f logs/bulk_atac_seq_*.out
```

---

## 📋 Expected Timeline and Results

### **Setup Time**: 
- Environment creation: 2-3 minutes
- Tool installation: 5-10 minutes
- BWA index creation: 10-20 minutes (first run only)

### **Processing Time per Sample**:
- Alignment: 1-3 hours
- Filtering and peak calling: 30 minutes - 1 hour
- **Total per sample**: 2-4 hours

### **Expected Results**:
```bash
# Check completion
grep "Processing pipeline completed" logs/bulk_atac_seq_*.out

# Verify output structure
ls -la bulk_atac_output/
# Should show:
# R26-Nestin-Ctrl-adult/
# R26-Nestin-Mut-adult/

# Check key files
ls -la bulk_atac_output/*/alignment/*.bam
ls -la bulk_atac_output/*/peaks/*.narrowPeak
ls -la bulk_atac_output/*/stats/*.txt
```

### **Quality Indicators**:
```bash
# Check alignment rates (should be >70%)
grep "mapped (" bulk_atac_output/*/stats/*_flagstat.txt

# Check peak counts (should be 20,000-100,000)
wc -l bulk_atac_output/*/peaks/*_peaks.narrowPeak
```

---

## 🛠️ Manual Commands (if needed)

### **Manually activate environment**:
```bash
conda activate atac_analysis
which bwa samtools macs2  # Verify tools
```

### **Manual tool installation** (if setup script fails):
```bash
# Create environment manually
conda create -n atac_analysis python=3.10 -y
conda activate atac_analysis
conda install -c bioconda -c conda-forge bwa samtools macs2 bowtie2 -y
```

### **Manual processing steps** (if batch job fails):
```bash
conda activate atac_analysis

# Process Control sample
cd bulk_atac_output
mkdir -p R26-Nestin-Ctrl-adult/alignment

# Align reads
bwa index /path/to/genome.fa  # if not done
bwa mem -t 16 /path/to/genome.fa \
    nestin/R26-Nestin-Ctrl-adult_R1_001.fastq.gz \
    nestin/R26-Nestin-Ctrl-adult_R3_001.fastq.gz \
    > R26-Nestin-Ctrl-adult/alignment/aligned.sam

# Convert and filter
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam -o sorted.bam
samtools view -b -f 2 -q 30 sorted.bam > filtered.bam

# Call peaks
macs2 callpeak -t filtered.bam -f BAMPE -g mm \
    -n R26-Nestin-Ctrl-adult --outdir peaks/
```

---

## 🚨 Troubleshooting

### **Environment creation fails**:
```bash
# Check conda installation
conda --version

# Update conda
conda update conda

# Try creating environment manually
conda create -n atac_analysis python=3.10 -y
```

### **Tool installation fails**:
```bash
# Add channels explicitly
conda config --add channels bioconda
conda config --add channels conda-forge

# Install tools individually
conda activate atac_analysis
conda install -c bioconda bwa -y
conda install -c bioconda samtools -y
conda install -c bioconda macs2 -y
```

### **Job fails to start**:
```bash
# Check if environment exists
conda env list | grep atac_analysis

# Check conda initialization
source ~/miniconda3/etc/profile.d/conda.sh
conda activate atac_analysis
```

### **Low alignment rates**:
- Check input files are correct (R1 and R3 genomic reads)
- Verify reference genome path
- Check for file corruption

---

## 🎯 Success Criteria

✅ **Environment created successfully** (`conda env list` shows `atac_analysis`)  
✅ **Tools available** (`which bwa samtools macs2` returns paths)  
✅ **Job completes** (logs show "Processing pipeline completed")  
✅ **BAM files generated** (>1GB each in alignment directories)  
✅ **Peaks called** (20,000+ peaks per sample)  
✅ **Good alignment rate** (>70% mapped reads)  

---

## 🎉 Next Steps After Success

1. **Load results into R/Python** for differential analysis
2. **Compare peak profiles** between Control vs Mutant
3. **Integrate with your RNA data** using standard bulk methods
4. **Generate publication figures** with tools like deepTools, IGV

---

**This approach completely bypasses the single-cell barcode issues and provides a robust solution for bulk ATAC-seq analysis of your samples.**