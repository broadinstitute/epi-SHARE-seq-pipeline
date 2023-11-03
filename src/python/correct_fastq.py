#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correct fastq
"""

import argparse
import xopen
from collections import deque


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


def create_barcode_dicts(barcode_list):
    """
    Creates dictionaries containing exact match and 1-mismatch sequences
    """
    barcode_dict_exact  = dict() # [seq: seq]
    barcode_dict_mismatch = dict() # [mismatch: seq]
    for barcode in barcode_list:
        barcode_dict_exact[barcode] = barcode # exact match
        for i, base in enumerate(barcode):
            for x in 'ACGTN':
                if base != x:
                    # add mismatch possibilities at pos i
                    mismatch = barcode[:i] + x + barcode[i + 1:]
                    barcode_dict_mismatch[mismatch] = barcode
                    
    return barcode_dict_exact, barcode_dict_mismatch

def check_putative_barcode(barcode_str, barcode_dict_exact, barcode_dict_mismatch, quality_str):
    """
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right
     shift
    """
    value = barcode_dict_exact.get(barcode_str[1:9]) # check exact location first
    quality = quality_str[1:9]
    if value is None:
        value = barcode_dict_mismatch.get(barcode_str[1:9])
        if value is None:
            value = barcode_dict_exact.get(barcode_str[:8]) # check 1bp shift left
            quality = quality_str[:8]
            if value is None:
                # check 1bp shift right
                # round 3 is shorter so add "N" for those
                if len(barcode_str) < 10: 
                    value = barcode_dict_exact.get(barcode_str[2:]+"N")
                    quality = quality_str[2:]+"F"
                else:
                    value = barcode_dict_exact.get(barcode_str[2:])
                    quality = quality_str[2:]
    return value, quality


def process_fastqs(input_read1_fastq_file,
                   input_read2_fastq_file,
                   output_read1_fastq_file,
                   output_read2_fastq_file,
                   output_barcode_fastq_file,
                   r1_exact_barcode_dict,
                   r1_mismatch_barcode_dict,
                   r2_exact_barcode_dict,
                   r2_mismatch_barcode_dict,
                   r3_exact_barcode_dict,
                   r3_mismatch_barcode_dict,
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
    # QC counters
    cellbarcode_match = cellbarcode_mismatch = cellbarcode_poly_g = read2_start_poly_g = 0

    read1_out_writer = xopen.xopen(output_read1_fastq_file, mode='w')
    read2_out_writer = xopen.xopen(output_read2_fastq_file, mode='w')
    barcode_out_writer = xopen.xopen(output_barcode_fastq_file, mode='w')

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
            r1_str, r2_str, r3_str = read_2_barcode_sequence[14:24], read_2_barcode_sequence[52:62], read_2_barcode_sequence[90:99]
            q1_str, q2_str, q3_str = read_2_barcode_quality[14:24], read_2_barcode_quality[52:62], read_2_barcode_quality[90:99]
            # get corrected barcodes
            r1 = r2 = r3 = None
            r1, q1 = check_putative_barcode(r1_str, r1_exact_barcode_dict, r1_mismatch_barcode_dict, q1_str)
            r2, q2 = check_putative_barcode(r2_str, r2_exact_barcode_dict, r2_mismatch_barcode_dict, q2_str)
            r3, q3 = check_putative_barcode(r3_str, r3_exact_barcode_dict, r3_mismatch_barcode_dict, q3_str)

            # check first ten base pairs of read 2 for homopolymer G
            if sequence2[:10] == "G"*10:
                read2_start_poly_g += 1

            # if corrected barcodes found, write to both read 1 and read 2 FASTQ files
            elif r1 and r2 and r3:
                cellbarcode_match += 1
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
                        corrected_quality2 = q1 + q2 + q3 + quality2[:-99]
                    else:
                        corrected_sequence2 = r1 + r2 + r3 + sequence2[:10]
                        corrected_quality2 = q1 + q2 + q3 + quality2[:10]
                    corrected_read2 = f"{corrected_header}\n{corrected_sequence2}\n+\n{corrected_quality2}\n"
                    buffer2.append(corrected_read2)
                    buffer_counter += 1
                elif sample_type == "ATAC":
                    # add corrected barcodes and PKR to header; remove any information after a space
                    corrected_header = name1.split(" ")[0] + "_" + ",".join(filter(None, [r1, r2, r3, pkr]))
                    # add corrected read 1 to buffer; use corrected header
                    corrected_read1 = f"{corrected_header}\n{sequence1}\n+\n{quality1}\n"
                    buffer1.append(corrected_read1)
                    # add corrected read 2 to buffer; use corrected header, remove 99bp barcode
                    corrected_read2 = f"{corrected_header}\n{sequence2[:-99]}\n+\n{quality2[:-99]}\n"
                    buffer2.append(corrected_read2)
                    #  add corrected barcode to buffer3 AKA fastq barcode
                    corrected_sequence3 = r1 + r2 + r3
                    corrected_quality3 = q1 + q2 + q3
                    corrected_read3 = f"{corrected_header}\n{corrected_sequence3}\n+\n{corrected_quality3}\n"
                    buffer3.append(corrected_read3)
                    buffer_counter += 1

                # write to corrected FASTQ files
                if buffer_counter == 10000000:
                    read1_out_writer.write("".join(buffer1))
                    buffer1.clear()
                    read2_out_writer.write("".join(buffer2))
                    buffer2.clear()
                    barcode_out_writer.write("".join(buffer3))
                    buffer3.clear()
                    buffer_counter = 0

            # check for homopolymer G in uncorrected barcode windows
            elif "G"*8 in r1_str and "G"*8 in r2_str and "G"*8 in r3_str:
                cellbarcode_poly_g += 1
            else:
                cellbarcode_mismatch += 1

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
        fields = ["library", "match", "mismatch", "poly_G_barcode", "poly_G_in_first_10bp"]
        f.write("\t".join(fields) + "\n")
        f.write("%s\t%s\t%s\t%s\t%s" % (prefix, cellbarcode_match, cellbarcode_mismatch, cellbarcode_poly_g, read2_start_poly_g))


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
    r1_exact_barcode_dict, r1_mismatch_barcode_dict = create_barcode_dicts(r1_barcodes)
    r2_exact_barcode_dict, r2_mismatch_barcode_dict = create_barcode_dicts(r2_barcodes)
    r3_exact_barcode_dict, r3_mismatch_barcode_dict = create_barcode_dicts(r3_barcodes)

    # write corrected FASTQs and QC stats
    process_fastqs(input_read1_fastq_file,
                   input_read2_fastq_file,
                   output_read1_fastq_file,
                   output_read2_fastq_file,
                   output_barcode_fastq_file,
                   r1_exact_barcode_dict,
                   r1_mismatch_barcode_dict,
                   r2_exact_barcode_dict,
                   r2_mismatch_barcode_dict,
                   r3_exact_barcode_dict,
                   r3_mismatch_barcode_dict,
                   sample_type,
                   pkr,
                   prefix,
                   paired_rna)


if __name__ == "__main__":
    main()
