# Spatial FFPE-GT-seq RNA Sequencing Data Processing

This part provides a workflow for generating RNA type-specific expression matrices from Spatial FFPE-GT-seq RNA sequencing data.

## Overview

The workflow consists of two major steps:

### Step 1: Generate Gene Expression Matrix using ASTRO

Spatial transcriptomics data are first processed using **ASTRO** to generate the complete RNA expression matrix.

ASTRO repository:

https://github.com/gersteinlab/ASTRO

The output of ASTRO serves as the input for downstream RNA type-specific analysis.

### Step 2: Split Expression Matrix by RNA Type

The script `split.R` separates the ASTRO-generated expression matrix into individual RNA type-specific matrices.

Examples of RNA types include:

* mRNA
* lncRNA
* snRNA
* snoRNA
* miRNA
* rRNA
* pseudogene
* other annotated transcript classes

Each RNA type is exported as an independent expression matrix for downstream analyses.

---


## Requirements

### Software

* R (≥ 4.3)
* ASTRO

### R Packages

```R
library(data.table)
library(dplyr)
```

Additional packages may be required depending on your implementation.

---
