# -*- coding: UTF-8 -*-
import re
import gzip
import os
import csv
import argparse
from multiprocessing import Pool
from natsort import natsorted
import pysam


def parse_args():
    """
    Parse command-line arguments for fragment counting.
    """
    parser = argparse.ArgumentParser(
        description="Count fragments per 5M base window from BAM/SAM files with UMI de-duplication"
    )
    parser.add_argument(
        '--directory', '-d', required=True,
        help='Directory containing sorted .bam files (with corresponding .bai) or .sam files'
    )
    parser.add_argument(
        '--fq2path', '-f', required=True,
        help='Path to the FASTQ file (unzipped .fq or gzipped .fq.gz)'
    )
    parser.add_argument(
        '--output_tsv', '-o', required=True,
        help='Path for output TSV file'
    )
    parser.add_argument(
        '--output_sum_csv', '-s', required=True,
        help='Path for output summary CSV file'
    )
    parser.add_argument(
        '--genome', '-g', choices=['human', 'mouse'], default='mouse',
        help='Select genome chromosome lengths: human or mouse'
    )
    return parser.parse_args()


def load_readsID_UMI(fq2path):
    """
    Load mapping of read IDs to UMI sequences from the FASTQ file.
    """
    readsID_UMI = {}
    pattern = re.compile(r'@(.*?)\s')
    open_func = gzip.open if fq2path.endswith('.gz') else open
    mode = 'rt' if fq2path.endswith('.gz') else 'r'
    with open_func(fq2path, mode) as fq2:
        for line_i, line in enumerate(fq2, start=1):
            if line_i % 4 == 1:
                m = pattern.search(line)
                if m:
                    thisKEY = m.group(1)
                else:
                    continue
            elif line_i % 4 == 2:
                # Extract UMI from positions 22-32 in read sequence line
                readsID_UMI[thisKEY] = line[22:32]
    return readsID_UMI


def get_chromosome_lengths(genome):
    """
    Return a dictionary of chromosome lengths for the specified genome.
    """
    if genome == 'human':
        return {
            'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
            'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
            'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
            'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
            'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
        }
    else:
        return {
            'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
            'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
            'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
            'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
            'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299, 'chrY': 91744698
        }


def count_fragments_per_5M_base(file_path, chromosome_lengths, window, readsID_UMI):
    """
    Count unique UMI-tagged fragments in 5M base windows across chromosomes.
    """
    # Initialize fragment counts per bin (set of UMIs)
    fragment_counts = {
        chrom: {i: set() for i in range((length // window) + 1)}
        for chrom, length in chromosome_lengths.items()
    }

    if file_path.endswith('.bam'):
        # Process BAM via pysam (requires existing .bai index)
        bam = pysam.AlignmentFile(file_path, "rb")
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            read_id = read.query_name
            chrom = read.reference_name
            pos = read.reference_start + 1
            umi = readsID_UMI.get(read_id)
            if not umi:
                continue
            bin_index = pos // window
            if chrom in fragment_counts and bin_index in fragment_counts[chrom]:
                fragment_counts[chrom][bin_index].add(umi)
        bam.close()
    else:
        # Fallback to SAM parsing
        with open(file_path, 'r') as fh:
            for line in fh:
                if line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                read_id, chrom, pos = parts[0], parts[2], parts[3]
                try:
                    pos = int(pos)
                except ValueError:
                    continue
                umi = readsID_UMI.get(read_id)
                if not umi:
                    continue
                bin_index = pos // window
                if chrom in fragment_counts and bin_index in fragment_counts[chrom]:
                    fragment_counts[chrom][bin_index].add(umi)

    # Flatten counts to list for output
    counts_list = []
    for chrom in natsorted(chromosome_lengths.keys()):
        num_bins = (chromosome_lengths[chrom] // window) + 1
        for i in range(num_bins):
            counts_list.append(len(fragment_counts[chrom][i]))
    return counts_list


def process_file(args_tuple):
    """
    Only process BAM files containing '+' in the filename; count fragments and return labeled array.
    """
    filename, directory, chrom_lengths, window, readsID_UMI = args_tuple
    if filename.endswith('.bam') and '+' in filename:
        filepath = os.path.join(directory, filename)

        # Check for .bai index
        bai1 = filepath + '.bai'
        bai2 = os.path.join(directory, filename[:-4] + '.bai')
        if not (os.path.exists(bai1) or os.path.exists(bai2)):
            print(f"[WARN] index for {filename} not found, skipping.")
            return None

        # Extract barcode pair (e.g., 'AACGTGAT+ACCACTGT') from filename
        base = os.path.splitext(filename)[0]
        m = re.search(r'([ACGT]{8}\+[ACGT]{8})', base)
        if not m:
            print(f"[WARN] cannot parse barcodes from {filename}, skipping.")
            return None
        part1, part2 = m.group(1).split('+')

        # Generate label using fixed mapping
        idx1 = 51 - int(cell_name_dic[part1])
        idx2 = int(cell_name_dic[part2])
        label = f"{idx1}x{idx2}"

        # Count fragments per 5M base bin
        result_array = count_fragments_per_5M_base(
            filepath, chrom_lengths, window, readsID_UMI
        )
        return [label] + result_array

    return None

# Fixed mapping from barcode to cell index
cell_name_dic = {
    'AACGTGAT': '1',   'AAACATCG': '2',   'ATGCCTAA': '3',   'AGTGGTCA': '4',   'ACCACTGT': '5',
    'ACATTGGC': '6',   'CAGATCTG': '7',   'CATCAAGT': '8',   'CGCTGATC': '9',   'ACAAGCTA': '10',
    'CTGTAGCC': '11',  'AGTACAAG': '12',  'AACAACCA': '13',  'AACCGAGA': '14',  'AACGCTTA': '15',
    'AAGACGGA': '16',  'AAGGTACA': '17',  'ACACAGAA': '18',  'ACAGCAGA': '19',  'ACCTCCAA': '20',
    'ACGCTCGA': '21',  'ACGTATCA': '22',  'ACTATGCA': '23',  'AGAGTCAA': '24',  'AGATCGCA': '25',
    'AGCAGGAA': '26',  'AGTCACTA': '27',  'ATCCTGTA': '28',  'ATTGAGGA': '29',  'CAACCACA': '30',
    'GACTAGTA': '31',  'CAATGGAA': '32',  'CACTTCGA': '33',  'CAGCGTTA': '34',  'CATACCAA': '35',
    'CCAGTTCA': '36',  'CCGAAGTA': '37',  'CCGTGAGA': '38',  'CCTCCTGA': '39',  'CGAACTTA': '40',
    'CGACTGGA': '41',  'CGCATACA': '42',  'CTCAATGA': '43',  'CTGAGCCA': '44',  'CTGGCATA': '45',
    'GAATCTGA': '46',  'GAAGACTA': '47',  'GAGCTGAA': '48',  'GATAGACA': '49',  'GCCACATA': '50'
}


def main():
    args = parse_args()
    directory = args.directory
    readsID_UMI = load_readsID_UMI(args.fq2path)
    chrom_lengths = get_chromosome_lengths(args.genome)
    window = 5000000

    # List only BAM files in directory
    filenames = [f for f in os.listdir(directory) if f.endswith('.bam')]
    task_args = [(fn, directory, chrom_lengths, window, readsID_UMI) for fn in filenames]
    with Pool() as pool:
        results = pool.map(process_file, task_args)

    filtered = [r for r in results if r]

    # Write detailed TSV output
    header = ['spot']
    for chrom in natsorted(chrom_lengths.keys()):
        num_bins = (chrom_lengths[chrom] // window) + 1
        header.extend([f"{chrom}_{i+1}" for i in range(num_bins)])
    with open(args.output_tsv, 'w', newline='') as out_tsv:
        writer = csv.writer(out_tsv, delimiter='\t')
        writer.writerow(header)
        writer.writerows(filtered)

    # Compute column sums for summary
    col_sums = [0] * (len(header) - 1)
    for row in filtered:
        for i, val in enumerate(row[1:], start=0):
            col_sums[i] += val

    # Write summary CSV
    with open(args.output_sum_csv, 'w', newline='') as out_sum:
        writer = csv.writer(out_sum)
        writer.writerow(['pos', 'fragments'])
        for idx, col in enumerate(header[1:], start=0):
            writer.writerow([col, col_sums[idx]])


if __name__ == '__main__':
    main()
