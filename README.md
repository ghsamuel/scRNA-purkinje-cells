# Purkinje Cell Diversity in Mouse Cerebellum Development

Single-cell RNA-seq analysis of Purkinje cell heterogeneity during embryonic mouse cerebellar development.

---

## Project Overview

This project analyzes single-cell RNA sequencing data from embryonic mouse cerebellum to investigate Purkinje cell diversity and developmental changes. The dataset (GSE256438) includes ~35,000 cells across two developmental timepoints, enabling identification of Purkinje cell subtypes and their molecular signatures.

**Dataset:** Khouri-Farah et al., *Nature Neuroscience* (2025)  
**DOI:** [10.1038/s41593-025-02042-w](https://doi.org/10.1038/s41593-025-02042-w)  
**GEO Accession:** GSE256438  
**BioProject:** PRJNA1079673

---

## Biological Question

How does Purkinje cell molecular heterogeneity emerge during cerebellar development between E16.5 and E18.5?

**Analysis Goals:**
1. Identify Purkinje cell subtypes using standard scRNA-seq clustering
2. Characterize developmental changes in cell-cell communication (CellChat)
3. **Extension:** Investigate alternative polyadenylation patterns in SMARCE1 across PC subtypes

---

## Dataset Details

| Sample | Age | Tech | Cells (est.) | Data Size |
|--------|-----|------|--------------|-----------|
| SRR28065512 | E16.5 rep1 | 10X Chromium | ~12K | ~18 GB |
| SRR28065511 | E16.5 rep2 | 10X Chromium | ~11K | ~17 GB |
| SRR28065510 | E18.5 | 10X Chromium | ~12K | ~19 GB |
| **Total** | - | - | **~35K** | **~54 GB** |

**Sequencing:** Illumina NextSeq 500, paired-end  
**Reference Genome:** mm10  
**Chemistry:** 10X Genomics v3

---

## Repository Structure
```
.
├── data/
│   └── fastq/              # Raw FASTQ files (not tracked in Git)
├── metadata/
│   ├── SraRunTable.csv     # SRA sample metadata
│   ├── sample_mapping.csv  # SRR → GSM → Stage mapping
│   └── subset.txt          # SRR accessions for download
├── scripts/
│   ├── 01_getData/
│   │   └── 02_sra.sh       # FASTQ download script
│   ├── 02_cellranger/      # Cell Ranger alignment
│   ├── 03_seurat/          # Seurat analysis scripts
│   └── 04_cellchat/        # CellChat analysis
├── results/                # Analysis outputs (not tracked)
├── logs/                   # SLURM logs
└── README.md
```

---

## Analysis Workflow

### 1. Data Download
```bash
# Download FASTQ files from SRA
sbatch scripts/01_getData/02_sra.sh
```

### 2. Quality Control
```bash
# FastQC on raw reads
# Cell Ranger web_summary.html
```

### 3. Alignment & Quantification
```bash
# Cell Ranger count per sample
# Output: feature-barcode matrices
```

### 4. Seurat Processing
- Quality filtering (nFeature, nCount, percent.mt)
- Sample integration (E16.5 rep1 + rep2)
- Dimensionality reduction (PCA, UMAP)
- Clustering and cell type annotation
- Differential expression analysis

### 5. CellChat Analysis
- Infer cell-cell communication networks
- Compare E16.5 vs E18.5 signaling patterns
- Identify developmental changes in ligand-receptor interactions

### 6. APA Extension (Optional)
- SMARCE1 isoform quantification across PC subtypes
- Integration with communication profiles

---

## Key Software

| Tool | Version | Purpose |
|------|---------|---------|
| Cell Ranger | 8.0+ | Alignment & quantification |
| R | 4.3+ | Analysis |
| Seurat | 5.0+ | scRNA-seq analysis |
| CellChat | 2.0+ | Cell communication inference |
| FastQC | 0.12+ | Quality control |
| MultiQC | 1.15+ | QC aggregation |

---

### Author

Glady Hazitha Samuel  
University of Connecticut
---

## References

Khouri-Farah, N., Guo, Q., Perry, T.A., Dussault, R., Li, J.Y.H. (2025). FOXP genes regulate Purkinje cell diversity and cerebellar morphogenesis. *Nature Neuroscience*. https://doi.org/10.1038/s41593-025-02042-w

---

## License

Analysis code: MIT License  
Data: GEO accession GSE256438 (published, public domain)
