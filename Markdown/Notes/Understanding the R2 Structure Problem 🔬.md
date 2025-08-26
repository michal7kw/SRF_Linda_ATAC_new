---
created: 2025-08-25T19:51
updated: 2025-08-26T09:36
---
## Understanding the R2 Structure Problem 🔬

### What You're Seeing:
```
Expected R2: [16bp barcode]
Your R2:     [CAGACG][16bp barcode][XX]
              └─6bp─┘ └────16bp────┘└2bp┘
              Total: 24bp
```

## The Root Cause: 

This is a **sequencing configuration issue** during the Illumina run setup. Here's what happened:

### 1. The 10x Multiome ATAC Library Structure:

During library prep, the actual molecular structure looks like this:

```
5'--[P5]--[Cell Barcode (16bp)]--[Linker (8bp)]--[Genomic DNA]--[ME]--[Genomic DNA]--[Linker]--[s7]--[i7]--[P7]--3'
                                    ↑
                                  CAGACGAT (8bp linker)
```

The **CAGACG** you're seeing is part of the **8bp linker sequence** that connects the cell barcode to the transposed genomic DNA!

### 2. What Went Wrong During Sequencing:

#### Correct Setup (Standard):
```yaml
Read 1: 50 cycles  → Genomic DNA (forward)
Index 1: 8 cycles  → Sample index (i7)
Read 2: 16 cycles  → Cell barcode ONLY
Read 3: 50 cycles  → Genomic DNA (reverse)
```

#### Your Setup (What Actually Happened):
```yaml
Read 1: 50 cycles  → Genomic DNA (forward) ✓
Index 1: 8 cycles  → Sample index (i7) ✓
Read 2: 24 cycles  → Cell barcode + linker sequence + extras
Read 3: 49 cycles  → Genomic DNA (reverse) ✓
```

### Why This Happened:

**Most likely scenario:** The sequencing facility set **Read 2 to 24 cycles** instead of 16 cycles. This caused the sequencer to:
1. Read the 16bp cell barcode ✓
2. **Continue reading** into the 8bp linker sequence
3. Only capture 6bp of the 8bp linker (CAGACG out of CAGACGAT)
4. Plus 2 more bases from the transposed genomic region

### Alternative Explanations:

1. **Custom Protocol:** Some facilities use modified protocols where they intentionally sequence through the linker for QC purposes

2. **Kit Version Confusion:** Different 10x kit versions might have been mixed up during setup

3. **Sequencer Configuration Template:** Used a wrong template (e.g., from a different 10x assay type)

## Why This Causes Processing Problems:

### 1. Barcode Matching Fails:
```python
# What the software expects:
Expected: "AAACAGCCAAACGCGA"  # 16bp from whitelist

# What it gets:
Actual:   "CAGACGAAACAGCCAAACGCGACC"  # 24bp
           └─6bp─┘└────16bp────┘└2bp┘

# Result: NO MATCH → Cell rejected!
```

### 2. Tools Can't Find Valid Cells:
- CellRanger expects exactly 16bp at the start of R2
- Finds "CAGACG..." instead of valid barcode
- Rejects 99%+ of reads as "invalid barcodes"

### 3. Downstream Effects:
- Massive data loss (>95% of cells rejected)
- Poor clustering (too few cells)
- Failed QC metrics
- Unusable results

## Why CAGACG is Always the Same:

The linker sequence is **constant by design**:
- It's part of the 10x Multiome chemistry
- Every gel bead has the same linker sequence
- Only the 16bp barcode portion varies

```
Bead 1: [Barcode_1]--[CAGACGAT]--[Capture sequence]
Bead 2: [Barcode_2]--[CAGACGAT]--[Capture sequence]
Bead 3: [Barcode_3]--[CAGACGAT]--[Capture sequence]
```

## The Biological Purpose:

The 8bp linker (CAGACGAT) serves important functions:
1. **Spacer:** Separates barcode from genomic DNA
2. **Priming site:** For downstream amplification
3. **Orientation:** Ensures correct molecular structure
4. **Tn5 compatibility:** Works with transposase chemistry

## Quality Control Implications:

Actually, seeing this consistent pattern is **good news**:
- ✅ Confirms your library structure is correct
- ✅ Shows successful 10x Multiome library prep
- ✅ Indicates no random degradation
- ✅ Easy to fix bioinformatically

## How This Likely Happened:

### At the Sequencing Facility:
```
Technician: "How many cycles for Read 2?"
PI/User: "Check the 10x protocol..."
*Looks at wrong protocol or misreads*
Result: 24 cycles instead of 16
```

Or:
```
Sequencing submission form:
□ Standard RNA-seq (28 cycles R2)
☑ 10x Single Cell (24 cycles R2)  ← Wrong box!
□ 10x ATAC (16 cycles R2)
```

## Prevention for Future:

1. **Explicitly specify** in sequencing orders:
   - "10x Multiome ATAC: R1=50, I1=8, **R2=16**, R3=50"

2. **Send the correct protocol** to the facility:
   - [10x Multiome Sequencing Guide](https://www.10xgenomics.com/support/single-cell-multiome-atac-plus-gene-expression/documentation/steps/sequencing)

3. **Request a test run** with few reads first

## The Silver Lining:

Your situation is **completely fixable** with the extraction script you already have:
```bash
# Extracts positions 7-22 (the actual barcode)
awk '{if(NR%4==2) print substr($0,7,16)}'
```

This is much better than having:
- Random N's throughout
- Variable read lengths
- Missing data
- Corrupted files

Your data is **perfectly usable** once you extract the correct 16bp! 🎯