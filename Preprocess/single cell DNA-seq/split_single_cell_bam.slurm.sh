#!/bin/bash
#SBATCH --job-name=split_bam
#SBATCH --output=logs/split_bam_%A_%a.out
#SBATCH --error=logs/split_bam_%A_%a.err
#SBATCH --partition=day
#SBATCH --array=1-7234%10  
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G


BAM="./possorted_bam.bam"
BARCODES="./barcodes for single cell.txt"
OUTDIR="./per_cell_bam_new"

mkdir -p $OUTDIR

BC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BARCODES)

echo "Processing barcode $BC"

samtools view -b -d CB:$BC $BAM \
    | samtools sort -@ 2 -o $OUTDIR/${BC}.sorted.bam

samtools index $OUTDIR/${BC}.sorted.bam

echo "Finished barcode $BC"
