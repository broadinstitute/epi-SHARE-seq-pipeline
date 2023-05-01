#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correct fastq
"""

import argparse
import dnaio

def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform barcode error correction on read 2 FASTQ file; write corrected barcodes into read names of both read 1 and read 2 FASTQ files; generate QC statistics file.")
    parser.add_argument("input_r1_fastq_file", help="Filename for uncorrected input read 1 FASTQ file")
    parser.add_argument("input_r2_fastq_file", help="Filename for uncorrected input read 2 FASTQ file")
    parser.add_argument("output_r1_fastq_file", help="Filename for corrected output read 1 FASTQ file")
    parser.add_argument("output_r2_fastq_file", help="Filename for corrected output read 2 FASTQ file")
    parser.add_argument("whitelist_file", help="Filename for whitelisted combinations of R1R2R3 barcodes, one per line")
    parser.add_argument("sample_type", choices=["ATAC", "RNA"], help="Sample modality")
    parser.add_argument("pkr", help="PKR name")
    parser.add_argument("prefix", help="Prefix for naming output QC txt file")
    
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

def create_barcode_dict(barcode_list):
    """
    Adds each barcode and its mismatch possibilities to the dictionary
    """
    barcode_dict  = dict()
    for barcode in barcode_list:
        barcode_dict[barcode] = barcode # exact match
        for i, base in enumerate(barcode):
            for x in 'ACGTN':
                if base != x:
                    # add mismatch possibilities at pos i
                    mismatch = barcode[:i] + x + barcode[i + 1:]
                    barcode_dict[mismatch] = barcode
                    
    return barcode_dict

def check_putative_barcode(barcode_str, barcode_dict, quality_str):
    """
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right shift
    """
    value = barcode_dict.get(barcode_str[1:9]) # check exact location first
    quality = quality_str[1:9]
    if value is None:
        value = barcode_dict.get(barcode_str[:8]) # check 1bp shift left
        quality = quality_str[:8]
        if value is None:
            # check 1bp shift right
            # round 3 is shorter so add "N" for those
            if len(barcode_str) < 10: 
                value = barcode_dict.get(barcode_str[2:]+"N")
                quality = quality_str[2:]+"F"
            else:
                value = barcode_dict.get(barcode_str[2:])
                quality = quality_str[2:]
                    
    return value, quality

def process_fastqs(input_r1_fastq_file, input_r2_fastq_file,
                  output_r1_fastq_file, output_r2_fastq_file,
                  r1_barcode_dict, r2_barcode_dict, r3_barcode_dict,
                  sample_type, pkr, prefix):
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
    
    # process FASTQs together
    with dnaio.open(input_r1_fastq_file, input_r2_fastq_file) as reader, \
         dnaio.open(output_r1_fastq_file, output_r2_fastq_file, mode="w") as writer:
             
        for read_1, read_2 in reader:
            # last 99bp of read 2 contains barcode sequences
            read_2_barcode_sequence = read_2.sequence[-99:]
            read_2_barcode_quality = read_2.qualities[-99:]
            # extract 10bp sequence containing R1 barcode, 10bp sequence containing R2 barcode, 
            # 9bp sequence containing R3 barcode, and corresponding quality strings
            r1_str, r2_str, r3_str = read_2_barcode_sequence[14:24], read_2_barcode_sequence[52:62], read_2_barcode_sequence[90:99]
            q1_str, q2_str, q3_str = read_2_barcode_quality[14:24], read_2_barcode_quality[52:62], read_2_barcode_quality[90:99]
            # get corrected barcodes
            r1 = r2 = r3 = None
            r1, q1 = check_putative_barcode(r1_str, r1_barcode_dict, q1_str)
            r2, q2 = check_putative_barcode(r2_str, r2_barcode_dict, q2_str)
            r3, q3 = check_putative_barcode(r3_str, r3_barcode_dict, q3_str)
            
            # check first ten base pairs of read 2 for homopolymer G
            if read_2.sequence[:10] == "G"*10:
                read2_start_poly_g += 1
                
            # if corrected barcodes found, write to both R1 and R2 FASTQ files
            elif r1 and r2 and r3:
                cellbarcode_match +=1
                # correct FASTQ reads
                if sample_type == "RNA":
                    # add corrected barcodes, PKR, and UMI to header
                    corrected_header = read_1.name + "_" + ",".join([r1, r2, r3, pkr]) + "_" + read_2.sequence[:10]                
                    # create SequenceRecord for read 1; use corrected header
                    corrected_read_1 = dnaio.SequenceRecord(corrected_header, read_1.sequence, read_1.qualities)                
                    # create SequenceRecord for read 2; use corrected header, read has format R1R2R3UMI
                    corrected_read_2_sequence = r1 + r2 + r3 + read_2.sequence[:10]
                    corrected_read_2_quality = q1 + q2 + q3 + read_2.qualities[:10]
                    corrected_read_2 = dnaio.SequenceRecord(corrected_header, corrected_read_2_sequence, corrected_read_2_quality)
                    
                elif sample_type == "ATAC":
                    # add corrected barcodes and PKR to header 
                    corrected_header = read_1.name + "_" + ",".join([r1, r2, r3, pkr])
                    # create SequenceRecord object for read 1: use corrected header
                    corrected_read_1 = dnaio.SequenceRecord(corrected_header, read_1.sequence, read_1.qualities)
                    # create SequenceRecord object for read 2: use corrected header, remove 99bp barcode
                    corrected_read_2 = dnaio.SequenceRecord(corrected_header, read_2.sequence[:-99], read_2.qualities[:-99])
                # write to corrected FASTQ files
                writer.write(corrected_read_1, corrected_read_2)
                
            # check for homopolymer G in uncorrected barcode windows
            elif "G"*8 in r1_str and "G"*8 in r2_str and "G"*8 in r3_str:
                cellbarcode_poly_g += 1
                
            else:
                cellbarcode_mismatch += 1
    
    # write QC stats
    with open(f"{prefix}_barcode_qc.txt", "w") as f:
        fields = ["library", "match", "mismatch", "poly_G_barcode", "poly_G_in_first_10bp"]
        f.write("\t".join(fields) + "\n")
        f.write("%s\t%s\t%s\t%s\t%s" % (prefix, cellbarcode_match, cellbarcode_mismatch, cellbarcode_poly_g, read2_start_poly_g))

def main():
    args = parse_arguments()
    input_r1_fastq_file = getattr(args, "input_r1_fastq_file")
    input_r2_fastq_file = getattr(args, "input_r2_fastq_file")
    output_r1_fastq_file = getattr(args, "output_r1_fastq_file")
    output_r2_fastq_file = getattr(args, "output_r2_fastq_file")
    whitelist_file = getattr(args, "whitelist_file")
    sample_type = getattr(args, "sample_type")
    pkr = getattr(args, "pkr")
    prefix = getattr(args, "prefix")
    
    # read whitelist, get lists of barcodes
    (r1_barcodes, r2_barcodes, r3_barcodes) = get_barcodes(whitelist_file)
    
    # create dictionaries for exact barcode matches and barcode mismatches
    r1_barcode_dict = create_barcode_dict(r1_barcodes)
    r2_barcode_dict = create_barcode_dict(r2_barcodes)
    r3_barcode_dict = create_barcode_dict(r3_barcodes)

    # write corrected FASTQs and QC stats
    process_fastqs(input_r1_fastq_file, input_r2_fastq_file,
                   output_r1_fastq_file, output_r2_fastq_file,
                   r1_barcode_dict, r2_barcode_dict, r3_barcode_dict,
                   sample_type, pkr, prefix)

if __name__ == "__main__":
    main()
