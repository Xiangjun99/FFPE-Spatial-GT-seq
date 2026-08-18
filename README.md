# FFPE Spatial-GT-seq

This repository contains the code associated with the manuscript **"Spatial genome-transcriptome co-mapping reveals tumor heterogeneity in FFPE tissues."**

## Overview

**FFPE Spatial-GT-seq** is a spatial multi-omics platform for joint profiling of genomic DNA, RNA, non-coding RNAs, and, in an extended format, immune-related proteins from the same formalin-fixed paraffin-embedded (FFPE) tissue section at single-cell resolution.

The method enables spatially resolved characterization of genomic alterations and transcriptional states within clinically archived tissues, providing a framework for investigating tumor heterogeneity, clonal evolution, transcriptional reprogramming, and non-coding RNA programs while preserving their spatial context.

<img width="3468" height="1535" alt="scheme_for_github" src="https://github.com/user-attachments/assets/29790f5b-ac03-4c9b-8416-07fe5a06f68e" />

## Repository structure

```text
FFPE-Spatial-GT-seq/
├── Preprocess/
│   ├── Bulk DNA seq/
│   ├── single cell DNA-seq/
│   ├── FFPE Spatial-GT-seq DNA/
│   ├── FFPE Spatial-GT-seq RNA/
│   └── FFPE Spatial-GT-seq ADT/
│
├── Visualization/
│
└── metrics comparisons/
```

### Preprocessing

The [`Preprocess`](Preprocess/) directory contains pipelines and scripts for processing sequencing data generated in this study, including:

1. **FFPE Spatial-GT-seq DNA**
   Processing of the DNA modality from FFPE Spatial-GT-seq experiments.

2. **FFPE Spatial-GT-seq RNA**
   Processing of the RNA modality from FFPE Spatial-GT-seq experiments.

3. **FFPE Spatial-GT-seq ADT**
   Processing and quantification of antibody-derived tag (ADT) sequencing data for spatial protein profiling.

4. **Single-cell DNA-seq**
   Processing of single-cell DNA sequencing data used in this study.

5. **Bulk DNA-seq**
   Processing of bulk DNA sequencing data used for genomic analyses and comparisons.

Detailed usage instructions, software requirements, and input/output descriptions are provided in the corresponding subdirectories.

## Visualization and downstream analysis

The [`Visualization`](Visualization/) directory contains scripts used for downstream analysis and visualization of processed multi-omics data, including:

* DNA analysis
* RNA analysis
* ADT analysis
* weighted nearest neighbor (WNN) analysis
* integration of copy-number variation (CNV) and RNA profiles
* gene co-expression analysis

## Metrics comparisons

The [`metrics comparisons`](metrics%20comparisons/) directory contains data and scripts used for evaluating sequencing and profiling performance and for comparison with previously published technologies.

## Data modalities

FFPE Spatial-GT-seq supports the following molecular modalities:

| Modality       | Information profiled                          |
| -------------- | --------------------------------------------- |
| DNA            | Genomic alterations and copy-number variation |
| RNA            | Spatial gene expression                       |
| Non-coding RNA | Spatial non-coding RNA programs               |
| ADT            | Immune-related protein abundance              |

Together, these measurements enable integrative analysis of genomic, transcriptomic, and protein-level heterogeneity within FFPE tissue sections.

## Software requirements

Software requirements differ among the individual preprocessing and downstream analysis workflows.

Please refer to the README files within each subdirectory for the corresponding:

* software dependencies and versions
* input data formats
* parameter settings
* preprocessing procedures
* output files
* usage instructions

## Contact

For questions regarding the code or analysis workflows, please open an issue in this repository.
