#!/usr/bin/env python3
import gzip
import os
import argparse
import subprocess
from collections import defaultdict
from multiprocessing import Pool
from functools import partial

# --------------------------
# 1. BWA 比对模块 (preprocess1)
# --------------------------
def run_bwa(fq1, ref, output_dir, threads):
    os.makedirs(output_dir, exist_ok=True)
    sam_path = os.path.join(output_dir, "aligned.sam")
    print(f"[1/4] Running BWA MEM alignment...\nOutput: {sam_path}")
    cmd = [
        "bwa", "mem", "-t", str(threads), ref, fq1
    ]
    with open(sam_path, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    print("[✔] Alignment finished.")
    return sam_path


# --------------------------
# 2. Barcode 提取与纠错 (preprocess2)
# --------------------------
def get_candidate_barcodes(barcode_file):
    with open(barcode_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def is_one_char_different(s1, s2):
    if len(s1) != len(s2):
        return False
    return sum(a != b for a, b in zip(s1, s2)) == 1

def process_fastq_chunk(lines, candidate_barcodes):
    cell_barcodes = defaultdict(list)
    for i in range(0, len(lines), 4):
        try:
            read_name = lines[i].strip()
            sequence = lines[i + 1].strip()

            barcode1 = sequence[76:84]
            barcode2 = sequence[38:46]
            barcode3 = sequence[0:8]

            corrected = []
            for bc in (barcode1, barcode2, barcode3):
                if bc in candidate_barcodes:
                    corrected.append(bc)
                else:
                    match = next((cb for cb in candidate_barcodes if is_one_char_different(bc, cb)), None)
                    corrected.append(match)

            if all(corrected):
                cell_bc = ''.join(corrected)
                clean_name = read_name.lstrip('@').split(' ')[0]
                cell_barcodes[cell_bc].append(clean_name)
        except IndexError:
            break
    return cell_barcodes

def extract_barcodes(fq2_file, candidate_barcodes, num_processes):
    cell_barcodes = defaultdict(list)
    lines = []
    with gzip.open(fq2_file, 'rt') as fq2:
        for line in fq2:
            lines.append(line)
            if len(lines) >= 400000:
                chunks = [lines[i:i+4000] for i in range(0, len(lines), 4000)]
                with Pool(num_processes) as pool:
                    results = pool.starmap(process_fastq_chunk, [(chunk, candidate_barcodes) for chunk in chunks])
                for res in results:
                    for bc, reads in res.items():
                        cell_barcodes[bc].extend(reads)
                lines = []
        if lines:
            chunks = [lines[i:i+4000] for i in range(0, len(lines), 4000)]
            with Pool(num_processes) as pool:
                results = pool.starmap(process_fastq_chunk, [(chunk, candidate_barcodes) for chunk in chunks])
            for res in results:
                for bc, reads in res.items():
                    cell_barcodes[bc].extend(reads)
    return cell_barcodes

def filter_barcodes_by_elbow(cell_barcodes, threshold=100):
    print(f"Barcode count before filtering: {len(cell_barcodes)}")
    filtered = {bc: reads for bc, reads in cell_barcodes.items() if len(reads) >= threshold}
    print(f"Barcode count after filtering: {len(filtered)}")
    return filtered

def build_read_to_barcode_mapping(cell_barcodes):
    return {read: bc for bc, reads in cell_barcodes.items() for read in reads}

def process_sam_chunk(chunk, read_to_bc):
    local_map = defaultdict(list)
    for line in chunk:
        if line.startswith('@'):
            continue
        read_name = line.split('\t')[0]
        if read_name in read_to_bc:
            local_map[read_to_bc[read_name]].append(line)
    return local_map

def split_sam_by_barcodes(sam_file, cell_barcodes, output_dir, num_processes):
    split_dir = os.path.join(output_dir, "split_sam")
    os.makedirs(split_dir, exist_ok=True)
    header_lines = []
    chunks = []
    chunk_size = 100000
    with open(sam_file, 'r') as sf:
        buf = []
        for line in sf:
            if line.startswith('@'):
                header_lines.append(line)
            else:
                buf.append(line)
                if len(buf) >= chunk_size:
                    chunks.append(buf)
                    buf = []
        if buf:
            chunks.append(buf)
    read_to_bc = build_read_to_barcode_mapping(cell_barcodes)
    with Pool(num_processes) as pool:
        results = pool.starmap(process_sam_chunk, [(c, read_to_bc) for c in chunks])
    final_map = defaultdict(list)
    for local in results:
        for bc, reads in local.items():
            final_map[bc].extend(reads)
    for bc, reads in final_map.items():
        out_file = os.path.join(split_dir, f"{bc}.sam")
        with open(out_file, 'w') as out:
            out.writelines(header_lines)
            out.writelines(reads)
    print("[✔] SAM files split by barcode.")
    return [os.path.join(split_dir, f) for f in os.listdir(split_dir) if f.endswith(".sam")]


# --------------------------
# 3. SAM → BAM 转换模块（优化版）
# --------------------------
def sam_to_bam(sam_files, threads, output_dir):
    bam_dir = os.path.join(output_dir, "bam")
    os.makedirs(bam_dir, exist_ok=True)
    print("[4/4] Converting SAM → sorted BAM and indexing...")
    
    for sam in sam_files:
        bc_name = os.path.basename(sam).replace(".sam", "")
        bam = os.path.join(bam_dir, f"{bc_name}.bam")
        
        # 直接 SAM → 排序 BAM
        subprocess.run([
            "samtools", "sort",
            "-@", str(threads),
            "-o", bam,
            sam
        ], check=True)
        
        # 建立索引
        subprocess.run(["samtools", "index", bam], check=True)
        
        # 删除原 SAM 文件
        os.remove(sam)
    
    print("[✔] All BAMs sorted and indexed, ready for CopyKit.")
    return bam_dir


# --------------------------
# 主函数入口
# --------------------------
def main():
    parser = argparse.ArgumentParser(description="Single-cell DNA preprocessing pipeline for CopyKit input.")
    parser.add_argument("--fq1", required=True, help="FASTQ file for alignment")
    parser.add_argument("--fq2", required=True, help="FASTQ file containing barcode information")
    parser.add_argument("--ref", required=True, help="Reference genome (FASTA)")
    parser.add_argument("--barcode_file", required=True, help="Candidate barcode list file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    fq1, fq2, ref = args.fq1, args.fq2, args.ref
    barcode_file, out_dir, threads = args.barcode_file, args.output, args.threads

    sam_path = run_bwa(fq1, ref, out_dir, threads)

    print("[2/4] Extracting and correcting cell barcodes...")
    candidates = get_candidate_barcodes(barcode_file)
    cell_bcs = extract_barcodes(fq2, candidates, threads)
    cell_bcs = filter_barcodes_by_elbow(cell_bcs)

    print("[3/4] Splitting SAM by barcodes...")
    sam_files = split_sam_by_barcodes(sam_path, cell_bcs, out_dir, threads)

    bam_dir = sam_to_bam(sam_files, threads, out_dir)
    print(f"\n✅ Pipeline completed successfully.\nBAM files ready in: {bam_dir}")


if __name__ == "__main__":
    main()

