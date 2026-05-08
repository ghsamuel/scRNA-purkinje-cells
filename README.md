# Purkinje Cell Diversity During Cerebellar Development

**Single-cell RNA-seq analysis pipeline for investigating Purkinje cell heterogeneity in embryonic mouse cerebellum**

[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41593--025--02042--w-blue)](https://doi.org/10.1038/s41593-025-02042-w)
[![GEO](https://img.shields.io/badge/GEO-GSE256438-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256438)
[![Docker](https://img.shields.io/badge/Docker-Containerized-2496ED?logo=docker)](Dockerfile)

---

## Overview

This repository contains a reproducible analysis pipeline for characterizing Purkinje cell diversity during mouse cerebellar development using single-cell RNA sequencing. The analysis integrates quality control, batch correction, cell type annotation, high-resolution subtyping, cell-cell communication inference, and alternative polyadenylation quantification.

**Dataset:** Khouri-Farah et al., *Nature Neuroscience* (2025)  
**Technology:** 10X Genomics Chromium single-cell 3' RNA-seq  
**Samples:** 19,723 high-quality cells across E16.5 (2 replicates) and E18.5  
**Reference Genome:** GRCm39 (mm39)

---

## Biological Question

How do Purkinje cells diversify during embryonic cerebellar development, and what cell-cell communication patterns coordinate this process?

**Key Findings:**
- Identified 15 Purkinje cell subtypes with distinct Foxp1/Foxp2/Foxp4 expression patterns
- Detected 15% reduction in cell-cell communication from E16.5 to E18.5
- Quantified 130 genes with significant 3'UTR length changes during development

---

## Dataset

| Sample | Stage | Replicate | Raw Reads | Cells (post-QC) | SRA Accession |
|--------|-------|-----------|-----------|-----------------|---------------|
| E16.5 rep1 | E16.5 | 1 | 286.2M | 6,409 | SRR28065512 |
| E16.5 rep2 | E16.5 | 2 | 415.3M | 6,692 | SRR28065511 |
| E18.5 | E18.5 | - | 378.4M | 6,622 | SRR28065510 |
| **Total** | - | - | **1,080M** | **19,723** | - |

**GEO Accession:** [GSE256438](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256438)  
**BioProject:** PRJNA1079673  
**Sequencing Platform:** Illumina NextSeq 500, paired-end  
**Chemistry:** 10X Genomics v3

---

## Repository Structure

```
scRNA-purkinje-cells/
├── scripts/
│   ├── 01_getData/           # SRA download scripts
│   ├── 02_qc/                # FastQC quality control
│   ├── 03_counts/            # Cell Ranger alignment
│   └── 04_R_analysis/        # Seurat analysis pipeline
│       ├── 01_data_loading_qc.R
│       ├── 02_integration_harmony.R
│       ├── 03_cell_type_annotation.R
│       ├── 04_purkinje_subtype_analysis.R
│       ├── 05_cellchat_analysis.R
│       └── 06_scalpel_apa_analysis.R
├── metadata/
│   ├── SraRunTable.csv       # Sample metadata
│   └── sample_mapping.csv    # SRR to stage mapping
├── Dockerfile                # Reproducible R environment
├── environment.yml           # Conda environment specification
├── purkinje_analysis_report.html  # Complete tutorial-style report
└── README.md
```

**Note:** Raw data (`data/`, `results/`) are not tracked in Git. See Data Availability below.

---

## Analysis Pipeline

### 1. Data Download & QC
- Download FASTQ files from SRA using `fasterq-dump --include-technical`
- Quality assessment with FastQC and MultiQC
- **Output:** Quality metrics, per-base sequence quality plots

### 2. Cell Ranger Alignment
- Build GRCm39 reference genome
- Align reads and generate feature-barcode matrices
- **Output:** `filtered_feature_bc_matrix/` per sample

### 3. Seurat Quality Control (Script 01)
- Load Cell Ranger outputs into Seurat
- Filter cells: 200-5,000 genes, >500 UMIs, <10% mitochondrial
- **Output:** 3 filtered Seurat objects (`.rds`)

### 4. Harmony Integration (Script 02)
- Merge samples and correct batch effects with Harmony
- Re-cluster on integrated embeddings
- **Output:** Integrated Seurat object, UMAP plots

### 5. Cell Type Annotation (Script 03)
- Find cluster markers with `FindAllMarkers()`
- Annotate 15 major cell types using canonical markers
- **Output:** Annotated Seurat object, marker gene tables

### 6. Purkinje Subtyping (Script 04)
- Extract Purkinje cells, re-cluster at high resolution (res = 1.0)
- Identify 15 subtypes based on Foxp1/Foxp2/Foxp4 expression
- **Output:** Purkinje subtype object, marker gene lists

### 7. CellChat Communication (Script 05)
- Infer cell-cell communication networks for E16.5 and E18.5
- Compare interaction strength and pathway usage
- **Output:** CellChat objects, network visualizations

### 8. Alternative Polyadenylation (Script 06)
- Run SCALPEL Nextflow pipeline on HPC
- Quantify 3'UTR length changes across development
- **Output:** Isoform counts, weighted 3'UTR lengths per gene

---

## Quick Start

### Option 1: Docker (Under Construction)

**Prerequisites:** Docker installed

```bash
# Clone repository
git clone https://github.com/ghsamuel/scRNA-purkinje-cells.git
cd scRNA-purkinje-cells

# Build Docker image
docker build -t purkinje-analysis .

# Download processed data (see Data Availability)
# Place RDS files in ~/data/

# Run analysis (scripts 02-05)
docker run -v ~/data:/data -v ~/results:/results purkinje-analysis \
  Rscript /app/scripts/02_integration_harmony.R

docker run -v ~/data:/data -v ~/results:/results purkinje-analysis \
  Rscript /app/scripts/03_cell_type_annotation.R
  
# Continue with scripts 04 and 05...
```

### Option 2: Conda Environment

```bash
# Create conda environment
conda env create -f environment.yml
conda activate purkinje-scrna

# Run scripts sequentially
Rscript scripts/04_R_analysis/01_data_loading_qc.R
Rscript scripts/04_R_analysis/02_integration_harmony.R
# ... etc
```

---

## Software Requirements

### Core Tools
- **R** 4.4.2
- **Cell Ranger** v9.0.1
- **Nextflow** 23.04+ (for SCALPEL)
- **Docker** or **Singularity/Apptainer** (for containerization)

### R Packages
- Seurat 5.4.0
- harmony 2.0.2
- CellChat 2.0+
- dplyr 1.2.1
- ggplot2 4.0.3
- patchwork 1.3.2

See `environment.yml` for complete list with exact versions.

---

## Data Availability

### Raw Data
- **FASTQ files:** [SRA PRJNA1079673](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1079673)

### Processed RDS Files (for Docker)
Filtered Seurat objects from Script 01 output will be available on Zenodo:  
🔗 **[DOI: 10.5281/zenodo.XXXXXXX]** *(to be added)*

---

## Results Summary

| Analysis | Key Finding |
|----------|-------------|
| **QC** | 19,723 / 24,434 cells retained (80.7%) |
| **Batch Correction** | Harmony successfully integrated E16.5 replicates |
| **Cell Types** | 15 major cell types identified |
| **Purkinje Subtypes** | 15 subtypes (vs. 11 in original paper) |
| **Communication** | 15% reduction E16.5 → E18.5 |
| **APA** | 130 genes with significant 3'UTR changes |

**Full results and figures:** See `purkinje_analysis_report.html`
---

## Troubleshooting

### Common Issues

**1. Missing Technical Reads**  
**Problem:** `fasterq-dump` discards barcodes/UMIs by default  
**Solution:** Use `--include-technical` flag

**2. Seurat v5 FindAllMarkers Error**  
**Problem:** Layers not joined  
**Solution:** Run `obj <- JoinLayers(obj)` before `FindAllMarkers()`




*Last updated: May 2026*
