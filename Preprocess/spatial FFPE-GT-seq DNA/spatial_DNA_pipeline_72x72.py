#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Full pipeline: BWA -> add CB/UB tags -> coordinate sort & index -> umi_tools dedup -> sort by CB -> per-CB BAM (sorted + indexed)
"""
import os
import argparse
import gzip
import pysam
import subprocess
from collections import defaultdict

parser = argparse.ArgumentParser(description="Single-cell DNA pipeline (one-script)")
parser.add_argument('-r1', '--fastq1', required=True, help='R1 FASTQ (gz ok)')
parser.add_argument('-r2', '--fastq2', required=True, help='R2 FASTQ (contains CB+UMI, gz ok)')
parser.add_argument('-ref', '--reference', required=True, help='Reference fasta for bwa mem')
parser.add_argument('-o', '--outdir', required=True, help='Output directory')
parser.add_argument('-t', '--threads', default=16, type=int, help='Number of threads for bwa/samtools')
parser.add_argument('--keep-intermediate', action='store_true', help='Do not delete intermediate files')
args = parser.parse_args()

R1 = args.fastq1
R2 = args.fastq2
REF = args.reference
OUTDIR = args.outdir
THREADS = args.threads
KEEP_INTER = args.keep_intermediate
os.makedirs(OUTDIR, exist_ok=True)

# file paths
aligned_bam       = os.path.join(OUTDIR, "aligned.bam")                      # bwa -> bam
tagged_bam        = os.path.join(OUTDIR, "aligned_tagged.bam")              # add CB/UB tags
sorted_bam        = os.path.join(OUTDIR, "aligned_tagged_sorted.bam")       # coord-sorted & indexed (input to dedup)
sorted_bam_idx    = sorted_bam + ".bai"
dedup_bam         = os.path.join(OUTDIR, "dedup.bam")                       # umi_tools output
dedup_sorted_cb   = os.path.join(OUTDIR, "dedup_sorted_byCB.bam")           # dedup -> sort by CB
per_spot_dir      = os.path.join(OUTDIR, "per_spot")                        # final per-CB BAMs
dedup_log         = os.path.join(OUTDIR, "dedup.log")

# barcode list (your set)
oneBarcodeArr = [
  'AACGTGAT','GCGAGTAA','ATGCCTAA','GCTAACGA','ACCACTGT','ACATTGGC','CAGATCTG','CATCAAGT',
  'CGCTGATC','ACAAGCTA','CTGTAGCC','AGTACAAG','AACAACCA','AACCGAGA','AACGCTTA','AAGACGGA',
  'AAGGTACA','ACACAGAA','ACAGCAGA','ACCTCCAA','ACGCTCGA','ACGTATCA','ACTATGCA','AGAGTCAA',
  'AGATCGCA','AGCAGGAA','AGTCACTA','ATCCTGTA','ATTGAGGA','CAACCACA','GACTAGTA','CAATGGAA',
  'CACTTCGA','CAGCGTTA','CATACCAA','CCAGTTCA','CCGAAGTA','CCGTGAGA','CCTCCTGA','CGAACTTA',
  'CGACTGGA','CGCATACA','CTCAATGA','CTGAGCCA','CTGGCATA','GAATCTGA','CAAGACTA','GAGCTGAA',
  'GATAGACA','GCCACATA','GCTCGGTA','GGAGAACA','GGTGCGAA','GTACGCAA','GTCGTAGA','GTCTGTCA',
  'GTGTTCTA','TAGGATGA','TATCAGCA','TCCGTCTA','TCTTCACA','TGAAGAGA','TGGAACAA','TGGCTTCA',
  'TGGTGGTA','TTCACGCA','AACTCACC','AAGAGATC','AAGGACAC','AATCCGTC','AATGTTGC','ACACGACC',
  'AAACATCG','AGTGGTCA'
]

def match_barcode(a, b, match_base_amount=7):
    return sum([1 for x, y in zip(a, b) if x == y]) >= match_base_amount

# -------------------------------
# 1) BWA MEM -> aligned BAM (direct to BAM)
# -------------------------------
print(f"[STEP 1] Running bwa mem -> BAM (threads={THREADS})")
# use process substitution, so keep shell=True and bash
cmd_bwa = f"bwa mem -t {THREADS} {REF} <(zcat {R1}) <(zcat {R2}) | samtools view -bS - > {aligned_bam}"
subprocess.run(cmd_bwa, shell=True, executable="/bin/bash", check=True)

# -------------------------------
# 2) Extract CB/UMI from R2 and create in-memory mapping
#    (simple exact/approx match according to provided barcodes)
# -------------------------------
print(f"[STEP 2] Extracting CB/UMI from R2 ({R2})")
readsID_barcode = {}
readsID_umi = {}
count = 0
with gzip.open(R2, "rt") as fh:
    line_i = 0
    rid = None
    for line in fh:
        line_i += 1
        if line_i % 4 == 1:
            rid = line.split()[0][1:]
        elif line_i % 4 == 2:
            seq = line.rstrip("\n")
            # positions used previously: UMI 22:32, B 32:40, A 70:78
            umi = seq[22:32]
            bc_b = seq[32:40]
            bc_a = seq[70:78] if len(seq) >= 78 else ""
            # fuzzy-match to known barcodes
            matched_b = None
            matched_a = None
            for bc in oneBarcodeArr:
                if matched_b is None and len(bc_b) == len(bc) and match_barcode(bc, bc_b):
                    matched_b = bc
                if matched_a is None and len(bc_a) == len(bc) and match_barcode(bc, bc_a):
                    matched_a = bc
                if matched_a and matched_b:
                    break
            if matched_a and matched_b:
                readsID_barcode[rid] = f"{matched_b}+{matched_a}"
                readsID_umi[rid] = umi
                count += 1
        if line_i % 1000000 == 0:
            print(f"  [INFO] processed {line_i} lines, mapped reads with CB+UMI so far: {count}")
print(f"[STEP 2] total reads with CB+UMI found: {count}")

# -------------------------------
# 3) Add CB/UB tags into BAM -> tagged_bam
# -------------------------------
print(f"[STEP 3] Adding CB/UB tags into BAM (writing {tagged_bam})")
in_bam = pysam.AlignmentFile(aligned_bam, "rb")
out_bam = pysam.AlignmentFile(tagged_bam, "wb", template=in_bam)

written = 0
for read in in_bam.fetch(until_eof=True):
    rid = read.query_name
    if rid in readsID_barcode:
        read.set_tag("CB", readsID_barcode[rid])
        read.set_tag("UB", readsID_umi[rid])
        out_bam.write(read)
        written += 1
    # optional: if you want all reads preserved even without CB, write them too
    # else: we only write reads that have CB mapping
    if written % 100000 == 0 and written > 0:
        print(f"  [INFO] wrote {written} tagged reads")
in_bam.close()
out_bam.close()
print(f"[STEP 3] finished writing {written} tagged reads to {tagged_bam}")

# -------------------------------
# 4) coordinate-sort tagged BAM and index (required for umi_tools)
# -------------------------------
print(f"[STEP 4] Sorting tagged BAM by coordinate -> {sorted_bam}")
subprocess.run(f"samtools sort -@ {THREADS} -o {sorted_bam} {tagged_bam}", shell=True, check=True)
subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)

if not KEEP_INTER:
    try:
        os.remove(aligned_bam)
    except Exception:
        pass
    try:
        os.remove(tagged_bam)
    except Exception:
        pass

# -------------------------------
# 5) umi_tools dedup (uses coordinate-sorted BAM)
# -------------------------------
print(f"[STEP 5] Running umi_tools dedup -> {dedup_bam}")
cmd_dedup = (
    f"umi_tools dedup "
    f"-I {sorted_bam} "
    f"-S {dedup_bam} "
    f"--extract-umi-method=tag "
    f"--umi-tag=UB "
    f"--cell-tag=CB "
    f"--log {dedup_log}"
)
subprocess.run(cmd_dedup, shell=True, check=True)

# -------------------------------
# 6) sort dedup bam by CB tag (prepare for splitting)
# -------------------------------
print(f"[STEP 6] Sorting dedup BAM by CB tag -> {dedup_sorted_cb}")
subprocess.run(f"samtools sort -@ {THREADS} -t CB -o {dedup_sorted_cb} {dedup_bam}", shell=True, check=True)

if not KEEP_INTER:
    try:
        os.remove(dedup_bam)
    except Exception:
        pass

# -------------------------------
# 7) Split per-CB: write unsorted per-CB BAMs, then sort & index each
# -------------------------------
print(f"[STEP 7] Splitting {dedup_sorted_cb} into per-CB BAMs in {per_spot_dir}")
os.makedirs(per_spot_dir, exist_ok=True)
bam_writers = defaultdict(lambda: None)

bam_in = pysam.AlignmentFile(dedup_sorted_cb, "rb")
for read in bam_in.fetch(until_eof=True):
    if read.has_tag("CB"):
        cb = read.get_tag("CB")
        if bam_writers[cb] is None:
            outpath = os.path.join(per_spot_dir, f"{cb}.unsorted.bam")
            bam_writers[cb] = pysam.AlignmentFile(outpath, "wb", template=bam_in)
        bam_writers[cb].write(read)
bam_in.close()

# close writers
for cb, w in bam_writers.items():
    if w:
        w.close()

# sort & index per-CB files
print("[STEP 7] Sorting and indexing each per-CB BAM (this may take time)...")
for fname in os.listdir(per_spot_dir):
    if fname.endswith(".unsorted.bam"):
        unsorted_path = os.path.join(per_spot_dir, fname)
        sorted_path = os.path.join(per_spot_dir, fname.replace(".unsorted.bam", ".bam"))
        subprocess.run(f"samtools sort -@ {THREADS} -o {sorted_path} {unsorted_path}", shell=True, check=True)
        subprocess.run(f"samtools index {sorted_path}", shell=True, check=True)
        if not KEEP_INTER:
            os.remove(unsorted_path)

print(f"[DONE] Pipeline finished. Per-CB BAMs available in {per_spot_dir}")
