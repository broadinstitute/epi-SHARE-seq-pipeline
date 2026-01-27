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
    parser.add_argument("whitelist_file", help="Filename for whitelisted combinations of R1R2R3 barcodes, one per line")
    parser.add_argument("sample_type", choices=["ATAC", "RNA"], help="Sample modality")
    parser.add_argument("prefix", help="Prefix for naming output QC txt file")
    parser.add_argument("pkr", nargs="?", help="PKR name")
    parser.add_argument("--barcode_fastq", help="Filename for separate barcode FASTQ file", default=None)
    
    return parser.parse_args()

def get_barcodes(whitelist_file):
    """
    Read barcode whitelist file, split into R1, R2, and R3 barcodes
    """
    r1_barcodes, r2_barcodes, r3_barcodes = set(), set(), set()
    with open(whitelist_file) as f:
        for line in f:
            r1_barcodes.add(line[:8])
            r2_barcodes.add(line[8:16])
            r3_barcodes.add(line[16:24])
    
    return r1_barcodes, r2_barcodes, r3_barcodes

def check_putative_barcode(barcode_str, barcode_set, quality_str):
    """
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right shift
    """

    # Helper function to find match (exact or 1-mismatch)
    def find_match(putative_barcode, barcode_set):
        # 1. Exact match
        if putative_barcode in barcode_set:
            return putative_barcode

        # 2. Mismatch search
        corrected_barcode = None
        # Generate neighbors of the *read* barcode
        for i, base in enumerate(putative_barcode):
            for x in 'ACGTN':
                if base != x:
                    neighbor = putative_barcode[:i] + x + putative_barcode[i + 1:]
                    if neighbor in barcode_set:
                        if corrected_barcode is not None and corrected_barcode != neighbor:
                            # Ambiguous match
                            return None
                        corrected_barcode = neighbor
        return corrected_barcode

    # Check exact location first
    value = find_match(barcode_str[1:9], barcode_set)
    quality = quality_str[1:9]
    if value is None:
        # Check 1bp shift left
        value = find_match(barcode_str[:8], barcode_set)
        quality = quality_str[:8]
        if value is None:
            # check 1bp shift right
            # round 3 is shorter so add "N" for those
            if len(barcode_str) < 10: 
                value = find_match(barcode_str[2:]+"N", barcode_set)
                quality = quality_str[2:]+"F"
            else:
                value = find_match(barcode_str[2:], barcode_set)
                quality = quality_str[2:]
                    
    return value, quality

def process_fastqs(input_read1_fastq_file, input_read2_fastq_file,
                  output_read1_fastq_file, output_read2_fastq_file,
                  r1_barcodes, r2_barcodes, r3_barcodes,
                  sample_type, pkr, prefix, barcode_fastq_file=None):
    """
    Takes in filenames for input and output FASTQ files, as well as
    dictionaries for R1, R2, R3 barcodes. 
    Corrects barcodes and writes corrected R1R2R3 sequence and corresponding quality
    string into output FASTQ files.
    For SHARE-seq, it expects barcodes in the last 99bp of Read 2, OR in the separate barcode_fastq_file if provided.
    Also produces txt file with barcode QC statistics; reports number of 
    exact barcode matches, non-exact barcode matches, non-matches, homopolymer G barcodes,
    homopolymer Gs in first 10bp of read 2 (UMI sequence for RNA, gDNA sequence for ATAC).
    """
    # QC counters
    cellbarcode_match = cellbarcode_mismatch = cellbarcode_poly_g = read2_start_poly_g = 0

    read1_out_writer = xopen.xopen(output_read1_fastq_file, mode = 'w')
    read2_out_writer = xopen.xopen(output_read2_fastq_file, mode ='w')

    buffer1 = deque()
    buffer2 = deque()
    buffer_counter = 0

    # Determine file openers
    read1_fh = xopen.xopen(input_read1_fastq_file, mode="r", threads=8)
    read2_fh = xopen.xopen(input_read2_fastq_file, mode="r", threads=8)
    barcode_fh = xopen.xopen(barcode_fastq_file, mode="r", threads=8) if barcode_fastq_file else None

    try:
        if barcode_fh:
             iterators = zip(read1_fh, read2_fh, barcode_fh)
        else:
             iterators = zip(read1_fh, read2_fh)

        for reads in iterators:
            readline1 = reads[0]
            readline2 = reads[1]
            readline_bc = reads[2] if barcode_fh else None

            name1 = readline1.strip()
            name2 = readline2.strip()

            # Advance iterators to get sequence
            readline1 = next(read1_fh)
            readline2 = next(read2_fh)
            readline_bc = next(barcode_fh) if barcode_fh else None

            sequence1 = readline1.strip()
            sequence2 = readline2.strip()
            sequence_bc = readline_bc.strip() if barcode_fh else None

            # Advance iterators to get +
            next(read1_fh)
            next(read2_fh)
            if barcode_fh: next(barcode_fh)

            # Advance iterators to get quality
            readline1 = next(read1_fh)
            readline2 = next(read2_fh)
            readline_bc = next(barcode_fh) if barcode_fh else None

            quality1 = readline1.strip()
            quality2 = readline2.strip()
            quality_bc = readline_bc.strip() if barcode_fh else None

            if barcode_fh:
                # Use barcode from separate file
                # If the separate file is shorter than 99bp, it might fail the slicing below
                # Assuming separate barcode file follows the same structure (last 99bp relevant or entire read)
                # But typically separate index reads are just the barcode.
                # HOWEVER, SHARE-seq structure described is 99bp containing 3 sub-barcodes.
                # If user provides a separate file, we assume it contains the barcode sequence.
                # If it's a raw 99bp read, we treat it same as `sequence2[-99:]`
                read_2_barcode_sequence = sequence_bc
                read_2_barcode_quality = quality_bc
                # If the barcode read is longer than 99bp, maybe take last 99?
                # Safest to assume the barcode READ is the barcode.
                if len(read_2_barcode_sequence) > 99:
                     read_2_barcode_sequence = read_2_barcode_sequence[-99:]
                     read_2_barcode_quality = read_2_barcode_quality[-99:]
            else:
                # last 99bp of read 2 contains barcode sequences
                read_2_barcode_sequence = sequence2[-99:]
                read_2_barcode_quality = quality2[-99:]

            # extract 10bp sequence containing R1 barcode, 10bp sequence containing R2 barcode, 
            # 9bp sequence containing R3 barcode, and corresponding quality strings
            if len(read_2_barcode_sequence) < 99:
                cellbarcode_mismatch += 1
                continue

            r1_str, r2_str, r3_str = read_2_barcode_sequence[14:24], read_2_barcode_sequence[52:62], read_2_barcode_sequence[90:99]
            q1_str, q2_str, q3_str = read_2_barcode_quality[14:24], read_2_barcode_quality[52:62], read_2_barcode_quality[90:99]

            # get corrected barcodes
            r1 = r2 = r3 = None
            r1, q1 = check_putative_barcode(r1_str, r1_barcodes, q1_str)
            r2, q2 = check_putative_barcode(r2_str, r2_barcodes, q2_str)
            r3, q3 = check_putative_barcode(r3_str, r3_barcodes, q3_str)
            
            # check first ten base pairs of read 2 for homopolymer G (UMI/gDNA)
            # This logic remains on sequence2 (the genomic read)
            if sequence2[:10] == "G"*10:
                read2_start_poly_g += 1
                
            # if corrected barcodes found, write to both read 1 and read 2 FASTQ files
            elif r1 and r2 and r3:
                cellbarcode_match +=1
                # correct FASTQ reads
                if sample_type == "RNA":
                    # add corrected barcodes, PKR, and UMI to header; remove any information after a space
                    corrected_header = name1.split(" ")[0] + "_" + ",".join(filter(None, [r1, r2, r3, pkr])) + "_" + sequence2[:10]                

                    # add corrected read 1 to buffer; use corrected header
                    corrected_read1 = f"{corrected_header}\n{sequence1}\n+\n{quality1}\n"
                    buffer1.append(corrected_read1)
                    # add corrected read 2 to buffer; use corrected header, read has format R1R2R3UMI
                    # If barcode file was provided, sequence2 is purely genomic (or whatever was passed as R2)
                    # The original logic appended `sequence2[:10]` (UMI).
                    # If barcode file is separate, we assume R2 is the cDNA read?
                    # "homopolymer Gs in first 10bp of read 2 (UMI sequence for RNA...)"
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

                    # add corrected read 2 to buffer; use corrected header
                    # Original logic: remove 99bp barcode from R2
                    if barcode_fh:
                        # If barcode is separate, R2 is likely the genomic read and shouldn't be trimmed
                        sequence2_out = sequence2
                        quality2_out = quality2
                    else:
                        sequence2_out = sequence2[:-99]
                        quality2_out = quality2[:-99]

                    corrected_read2 = f"{corrected_header}\n{sequence2_out}\n+\n{quality2_out}\n"
                    buffer2.append(corrected_read2)
                    buffer_counter += 1

                # write to corrected FASTQ files
                if buffer_counter == 10000000:
                    read1_out_writer.write("".join(buffer1))
                    buffer1.clear()
                    read2_out_writer.write("".join(buffer2))
                    buffer2.clear()
                    buffer_counter = 0
                
            # check for homopolymer G in uncorrected barcode windows
            elif "G"*8 in r1_str and "G"*8 in r2_str and "G"*8 in r3_str:
                cellbarcode_poly_g += 1
                
            else:
                cellbarcode_mismatch += 1

    finally:
        read1_fh.close()
        read2_fh.close()
        if barcode_fh:
            barcode_fh.close()

    if buffer_counter > 0:
        read1_out_writer.write("".join(buffer1))
        buffer1.clear()
        read2_out_writer.write("".join(buffer2))
        buffer2.clear()
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
    whitelist_file = getattr(args, "whitelist_file")
    sample_type = getattr(args, "sample_type")
    prefix = getattr(args, "prefix")
    pkr = getattr(args, "pkr")
    barcode_fastq_file = getattr(args, "barcode_fastq")
    
    # read whitelist, get lists of barcodes
    (r1_barcodes, r2_barcodes, r3_barcodes) = get_barcodes(whitelist_file)

    # write corrected FASTQs and QC stats
    process_fastqs(input_read1_fastq_file, input_read2_fastq_file,
                   output_read1_fastq_file, output_read2_fastq_file,
                   r1_barcodes, r2_barcodes, r3_barcodes,
                   sample_type, pkr, prefix, barcode_fastq_file)

if __name__ == "__main__":
    main()
