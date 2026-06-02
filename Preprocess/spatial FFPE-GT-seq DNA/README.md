# Spatial-GT-seq DNA Sequencing Data Processing

This part of the pipeline processes spatially barcoded bulk DNA sequencing data and outputs barcode-specific BAM files and fragment count matrices.

## Step 1: Preprocessing (`preprocess.py`)

This script performs:

* Read alignment using BWA MEM
* SAM to BAM conversion, sorting, and indexing using Samtools
* Duplicate marking and removal using Sambamba
* Splitting deduplicated SAM files by barcode (extracted from FASTQ read 2)
* Per-barcode SAM to BAM conversion and autosomal read count filtering

#### Usage

```bash
python preprocess.py \
  --fq1 path/to/sample_R1.fq.gz \
  --fq2 path/to/sample_R2.fq.gz \
  --ref path/to/reference.fa \
  --output path/to/output_directory \
  --threshold 5
```

The final filtered BAM files will be located in:

```
<output_directory>/bamFile/
```

### Step 2: Fragment Counting (`getFragmentsInfo.py`)

This script counts unique UMI-tagged fragments for each BAM file in 5 million base pair bins across all chromosomes.

#### Usage

```bash
python getFragmentsInfo.py \
  --directory path/to/bamFile \
  --fq2path path/to/sample_R2.fq.gz \
  --output_tsv fragments_matrix.tsv \
  --output_sum_csv genome_bin_summary.csv \
  --genome mouse (or human)
```

Output:

* `fragments_matrix.tsv`: Each row represents a barcode spot; each column is a genomic bin.
* `genome_bin_summary.csv`: Total fragment count per genomic bin.

---


## Output Explanation

* **BAM files (spatial)**: Ready for downstream CNV calling or genome-wide coverage analysis.
* **Fragment matrix**: Useful for clustering, dimensionality reduction, and signal heatmaps.

---

## Performance Optimization

* Scripts utilize Python's `multiprocessing.Pool` to maximize CPU utilization.
* Ensure adequate system memory and thread capacity when processing large datasets.

---

