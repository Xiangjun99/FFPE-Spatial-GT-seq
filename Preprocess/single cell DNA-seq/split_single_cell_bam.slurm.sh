#!/bin/bash
#SBATCH --job-name=split_bam
#SBATCH --output=logs/split_bam_%A_%a.out
#SBATCH --error=logs/split_bam_%A_%a.err
#SBATCH --partition=scavenge
#SBATCH --array=1-7234%10   # 这里 1-100 是 barcode 数量，你可以改，%10 表示同时运行 10 个任务
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G


# 定义输入输出
BAM="/gpfs/ycga/work/liu_yang/xd97/@@2026_coprofile_revised/12_single_cell_ATAC/SCKPN11/outs/possorted_bam.bam"
BARCODES="/gpfs/ycga/work/liu_yang/xd97/@@2026_coprofile_revised/12_single_cell_ATAC/split/valid_barcodes_no_prefix1.txt"
OUTDIR="/gpfs/ycga/work/liu_yang/xd97/@@2026_coprofile_revised/12_single_cell_ATAC/split/per_cell_bam_new"

mkdir -p $OUTDIR

# 获取当前任务对应的 barcode
BC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BARCODES)

echo "Processing barcode $BC"

# 拆分 BAM 并排序
samtools view -b -d CB:$BC $BAM \
    | samtools sort -@ 2 -o $OUTDIR/${BC}.sorted.bam

# 生成索引
samtools index $OUTDIR/${BC}.sorted.bam

echo "Finished barcode $BC"
