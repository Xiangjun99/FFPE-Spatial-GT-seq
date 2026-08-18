# Bulk DNA-seq Preprocessing

This directory contains a Snakemake workflow for preprocessing paired-end bulk DNA sequencing data.

The workflow performs:

1. Paired-end read alignment using Bowtie2
2. BAM conversion and sorting using SAMtools
3. Duplicate marking using Sambamba

The final output consists of duplicate-marked BAM files.

## Workflow

The workflow is implemented in:
bowtie2_paired.smk

---

## Input files

Paired-end FASTQ files should be gzip-compressed and follow the naming convention:

{sample}_R1.fastq.gz

{sample}_R2.fastq.gz

---

## Software versions

The workflow was run using the following software versions:

| Software  | Version | Purpose                               |
| --------- | ------: | ------------------------------------- |
| Snakemake |  8.29.2 | Workflow management                   |
| Bowtie2   |   2.5.4 | Paired-end read alignment             |
| SAMtools  |    1.21 | BAM conversion, sorting, and indexing |
| Sambamba  |   1.0.1 | Duplicate marking                     |

---

## Reference genomes

Bowtie2 genome indexes are required before running the workflow.

Two reference genome indexes were used depending on the species.

---

### Human reference genome

Species:

Homo sapiens

Genome assembly:

GRCh38

Bowtie2 index prefix:

GRCh38_noalt_decoy_as

---

### Mouse reference genome

Species:

Mus musculus

Genome assembly:
mm10

Bowtie2 index prefix:

mm10

---

## Output files

Intermediate alignment files are generated in:

mapped/

Sorted BAM files are generated in:

sort/

Final duplicate-marked BAM files are generated in:

marked/

For example:

marked/Sample01.bam

marked/Sample02.bam

---

## Log files

Bowtie2 alignment logs are written to:

logs/bowtie2/

For example:

logs/bowtie2/Sample01.log

These logs contain Bowtie2 alignment statistics and are useful for checking mapping performance and troubleshooting failed runs.

---
