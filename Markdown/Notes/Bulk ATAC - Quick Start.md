---
created: 2025-08-26T08:16
updated: 2025-08-26T09:13
---
## 🛠️ Tool Installation (Fixed for Python 3.12)

### Method 1: Step-by-step Installation (Recommended)
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data

# Install BWA and samtools via conda
conda install -c bioconda bwa samtools -y

# Install MACS2 via pip (avoids conda conflicts)
pip install MACS2

# Verify installations
which bwa samtools macs2
```

### Method 2: Use Alternative Installation Script
```bash
bash install_tools_alternative.sh
```

### Method 3: Manual Installation (if above fails)
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

**The key point is that BWA + samtools are sufficient to generate aligned BAM files, and peak calling can be done separately if needed.**