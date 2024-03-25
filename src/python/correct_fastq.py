#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correct fastq
"""

import argparse
import xopen
from collections import deque
from utils import check_putative_barcode, create_barcode_dicts


def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform barcode error correction on read 2 FASTQ file; write corrected barcodes into read names of both read 1 and read 2 FASTQ files; generate QC statistics file.")
    parser.add_argument("input_read1_fastq_file", help="Filename for uncorrected input read 1 FASTQ file")
    parser.add_argument("input_read2_fastq_file", help="Filename for uncorrected input read 2 FASTQ file")
    parser.add_argument("output_read1_fastq_file", help="Filename for corrected output read 1 FASTQ file")
    parser.add_argument("output_read2_fastq_file", help="Filename for corrected output read 2 FASTQ file")
    parser.add_argument("output_barcode_fastq_file", help="Filename for corrected output barcode FASTQ file")
    parser.add_argument("whitelist_file", help="Filename for whitelisted combinations of R1R2R3 barcodes, one per line")
    parser.add_argument("sample_type", choices=["ATAC", "RNA"], help="Sample modality")
    parser.add_argument("prefix", help="Prefix for naming output QC txt file")
    parser.add_argument("pkr", nargs="?", help="PKR name")
    parser.add_argument("--paired_rna", help="Flag for outputting the full RNA read 2 for paired-end alignment, rather than just CB+UMI", action="store_true")

    return parser.parse_args()


def get_barcodes(whitelist_file):
    """
    Read barcode whitelist file, split into R1, R2, and R3 barcodes
    """
    r1_barcodes = r2_barcodes = r3_barcodes = set()
    with open(whitelist_file) as f:
        for line in f:
            r1_barcodes.add(line[:8])
            r2_barcodes.add(line[8:16])
            r3_barcodes.add(line[16:24])

    return list(r1_barcodes), list(r2_barcodes), list(r3_barcodes)


def make_correction_counts_dict(sample_type):
    """
    Initialize a dictionary for counting cell barcode correction combinations.

    This function initializes a dictionary to count the different combinations of cell barcode corrections.
    Accepted combinations are those which do not contain shifts.
    The dictionary also contains a key for counting non-matches and, if sample_type is "RNA", a key for counting polyG UMIs.

    Args:
        sample_type (str): The type of sample. If "RNA", an additional key for counting polyG UMIs will be added.

    Returns:
        dict: A dictionary with keys representing different correction combinations and their counts.
            The dictionary also contains a key for counting non-matches and, if sample_type is "RNA", a key for counting polyG UMIs.
    """
    cellbarcode_correction_counts = {}
    correction_types = ["E", "M"]  # exact, mismatch, right shift, left shift
    correction_combos = [c1+c2+c3 for c1 in correction_types for c2 in correction_types for c3 in correction_types]

    for combo in correction_combos:
        cellbarcode_correction_counts[combo] = 0

    cellbarcode_correction_counts["nonmatch"] = 0
    if sample_type == "RNA":
        cellbarcode_correction_counts["polyG_umi"] = 0

    return cellbarcode_correction_counts


def process_fastqs(input_read1_fastq_file,
                   input_read2_fastq_file,
                   output_read1_fastq_file,
                   output_read2_fastq_file,
                   output_barcode_fastq_file,
                   r1_barcode_exact_dict,
                   r1_barcode_mismatch_dict,
                   r2_barcode_exact_dict,
                   r2_barcode_mismatch_dict,
                   r3_barcode_exact_dict,
                   r3_barcode_mismatch_dict,
                   sample_type,
                   pkr,
                   prefix,
                   paired_rna):
    """
    Takes in filenames for input and output FASTQ files, as well as
    dictionaries for R1, R2, R3 barcodes.
    Corrects barcodes and writes corrected R1R2R3 sequence and corresponding quality
    string into output FASTQ files.
    Also produces txt file with barcode QC statistics; reports number of
    exact barcode matches, non-exact barcode matches, non-matches, homopolymer G barcodes,
    homopolymer Gs in first 10bp of read 2 (UMI sequence for RNA, gDNA sequence for ATAC).
    """
    # initialize dictionary with cell barcode correction types
    cellbarcode_correction_counts = make_correction_counts_dict(sample_type)

    read1_out_writer = xopen.xopen(output_read1_fastq_file, mode="w")
    read2_out_writer = xopen.xopen(output_read2_fastq_file, mode="w")
    barcode_out_writer = xopen.xopen(output_barcode_fastq_file, mode="w")

    buffer1 = deque()
    buffer2 = deque()
    buffer3 = deque()
    buffer_counter = 0

    # process FASTQs together
    with xopen.xopen(input_read1_fastq_file, mode="r", threads=8) as read1_fh, xopen.xopen(input_read2_fastq_file, mode="r", threads=8) as read2_fh:
        for readline1, readline2 in zip(read1_fh, read2_fh):
            name1 = readline1.strip()
            name2 = readline2.strip()

            readline1 = next(read1_fh)
            readline2 = next(read2_fh)

            sequence1 = readline1.strip()
            sequence2 = readline2.strip()

            next(read1_fh)
            next(read2_fh)

            readline1 = next(read1_fh)
            readline2 = next(read2_fh)

            quality1 = readline1.strip()
            quality2 = readline2.strip()

            # last 99bp of read 2 contains barcode sequences
            read_2_barcode_sequence = sequence2[-99:]
            read_2_barcode_quality = quality2[-99:]
            # extract 10bp sequence containing R1 barcode, 10bp sequence containing R2 barcode,
            # 9bp sequence containing R3 barcode, and corresponding quality strings
            r1_str, r2_str, r3_str = read_2_barcode_sequence[15:23], read_2_barcode_sequence[53:61], read_2_barcode_sequence[91:99]
            q1_str, q2_str, q3_str = read_2_barcode_quality[15:23], read_2_barcode_quality[53:61], read_2_barcode_quality[91:99]
            # get corrected barcodes
            r1 = r2 = r3 = None
            r1, c1 = check_putative_barcode(r1_str, r1_barcode_exact_dict, r1_barcode_mismatch_dict)
            r2, c2 = check_putative_barcode(r2_str, r2_barcode_exact_dict, r2_barcode_mismatch_dict)
            r3, c3 = check_putative_barcode(r3_str, r3_barcode_exact_dict, r3_barcode_mismatch_dict)

            # if corrected barcodes found and correction type is valid,
            # write to both read 1 and read 2 FASTQ files
            if r1 and r2 and r3:
                # correct FASTQ reads
                if sample_type == "RNA":
                    # add corrected barcodes, PKR, and UMI to header; remove any information after a space
                    corrected_header = name1.split(" ")[0] + "_" + ",".join(filter(None, [r1, r2, r3, pkr])) + "_" + sequence2[:10]

                    # add corrected read 1 to buffer; use corrected header
                    corrected_read1 = f"{corrected_header}\n{sequence1}\n+\n{quality1}\n"
                    buffer1.append(corrected_read1)
                    # add corrected read 2 to buffer; use corrected header, read has format R1R2R3UMI if not paired-end RNA,
                    # or R1R2R3UMIread2 if paired-end RNA
                    if paired_rna:
                        corrected_sequence2 = r1 + r2 + r3 + sequence2[:-99]
                        corrected_quality2 = q1_str + q2_str + q3_str + quality2[:-99]
                    else:
                        corrected_sequence2 = r1 + r2 + r3 + sequence2[:10]
                        corrected_quality2 = q1_str + q2_str + q3_str + quality2[:10]
                    corrected_read2 = f"{corrected_header}\n{corrected_sequence2}\n+\n{corrected_quality2}\n"
                    buffer2.append(corrected_read2)
                    buffer_counter += 1

                    # check for polyG UMI
                    if sequence2[:10] == "G"*10:
                        cellbarcode_correction_counts["polyG_umi"] += 1

                elif sample_type == "ATAC":
                    # add corrected barcodes and PKR to header; remove any information after a space
                    corrected_header = name1.split(" ")[0] + "_" + ",".join(filter(None, [r1, r2, r3, pkr]))
                    # add corrected read 1 to buffer; use corrected header
                    corrected_read1 = f"{corrected_header}\n{sequence1}\n+\n{quality1}\n"
                    buffer1.append(corrected_read1)
                    # add corrected read 2 to buffer; use corrected header, remove 99bp barcode
                    corrected_read2 = f"{corrected_header}\n{sequence2[:-99]}\n+\n{quality2[:-99]}\n"
                    buffer2.append(corrected_read2)
                    # add corrected barcode to buffer3 AKA fastq barcode
                    corrected_sequence3 = r1 + r2 + r3
                    corrected_quality3 = q1_str + q2_str + q3_str
                    corrected_read3 = f"{corrected_header}\n{corrected_sequence3}\n+\n{corrected_quality3}\n"
                    buffer3.append(corrected_read3)
                    buffer_counter += 1

                # add to correction counts dictionary
                cellbarcode_correction_counts[c1+c2+c3] += 1

                # write to corrected FASTQ files
                if buffer_counter == 10000000:
                    read1_out_writer.write("".join(buffer1))
                    buffer1.clear()
                    read2_out_writer.write("".join(buffer2))
                    buffer2.clear()
                    barcode_out_writer.write("".join(buffer3))
                    buffer3.clear()
                    buffer_counter = 0
            else:
                cellbarcode_correction_counts["nonmatch"] += 1

    if buffer_counter > 0:
        read1_out_writer.write("".join(buffer1))
        buffer1.clear()
        read2_out_writer.write("".join(buffer2))
        buffer2.clear()
        barcode_out_writer.write("".join(buffer3))
        buffer3.clear()
        buffer_counter = 0

    # write QC stats
    with open(f"{prefix}_barcode_qc.txt", "w") as f:
        f.write(f"\t{prefix}\n")
        for k, v in cellbarcode_correction_counts.items():
            f.write(f"{k}\t{v}\n")


def main():
    args = parse_arguments()
    input_read1_fastq_file = getattr(args, "input_read1_fastq_file")
    input_read2_fastq_file = getattr(args, "input_read2_fastq_file")
    output_read1_fastq_file = getattr(args, "output_read1_fastq_file")
    output_read2_fastq_file = getattr(args, "output_read2_fastq_file")
    output_barcode_fastq_file = getattr(args, "output_barcode_fastq_file")
    whitelist_file = getattr(args, "whitelist_file")
    sample_type = getattr(args, "sample_type")
    prefix = getattr(args, "prefix")
    pkr = getattr(args, "pkr")
    paired_rna = getattr(args, "paired_rna")

    # read whitelist, get lists of barcodes
    (r1_barcodes, r2_barcodes, r3_barcodes) = get_barcodes(whitelist_file)

    # create dictionaries for exact barcode matches and barcode mismatches
    r1_barcode_exact_dict, r1_barcode_mismatch_dict = create_barcode_dicts(r1_barcodes)
    r2_barcode_exact_dict, r2_barcode_mismatch_dict = create_barcode_dicts(r2_barcodes)
    r3_barcode_exact_dict, r3_barcode_mismatch_dict = create_barcode_dicts(r3_barcodes)

    # write corrected FASTQs and QC stats
    process_fastqs(input_read1_fastq_file,
                   input_read2_fastq_file,
                   output_read1_fastq_file,
                   output_read2_fastq_file,
                   output_barcode_fastq_file,
                   r1_barcode_exact_dict,
                   r1_barcode_mismatch_dict,
                   r2_barcode_exact_dict,
                   r2_barcode_mismatch_dict,
                   r3_barcode_exact_dict,
                   r3_barcode_mismatch_dict,
                   sample_type,
                   pkr,
                   prefix,
                   paired_rna)


if __name__ == "__main__":
    main()
