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
    parser.add_argument("r1_barcode_file", help="File containing R1 barcodes, one line")
    parser.add_argument("r2_barcode_file", help="File containing R2 barcodes, one line")
    parser.add_argument("r3_barcode_file", help="File containing R3 barcodes, one line")
    parser.add_argument("sample_type", choices=["ATAC", "RNA"], help="Sample modality")
    parser.add_argument("pkr", help="PKR name")
    parser.add_argument("prefix", help="Prefix for naming output QC txt file")
    
    return parser.parse_args()

def create_barcode_dicts(barcode_list):
    """
    Adds each barcode and its mismatch possibilities to the dictionary
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
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right shift
    """
    exact = True
    value = barcode_dict_exact.get(barcode_str[1:9]) # check exact location first
    quality = quality_str[1:9]
    if value is None:
        exact = False
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
                    
    return value, quality, exact

def process_fastqs(input_r1_fastq_file, input_r2_fastq_file,
                  output_r1_fastq_file, output_r2_fastq_file,
                  r1_barcode_dict_exact, r1_barcode_dict_mismatch, 
                  r2_barcode_dict_exact, r2_barcode_dict_mismatch, 
                  r3_barcode_dict_exact, r3_barcode_dict_mismatch,
                  sample_type, pkr, prefix):
    """
    Takes in filenames for input and output FASTQ files, as well as exact and mismatch
    dictionaries for R1, R2, R3 barcodes. 
    Corrects barcodes and writes corrected R1R2R3 sequence and corresponding quality
    string into output FASTQ files.
    Also produces txt file with barcode QC statistics; reports number of 
    exact barcode matches, non-exact barcode matches, non-matches, homopolymer G barcodes,
    homopolymer Gs in first 10bp of read 2 (UMI sequence for RNA, gDNA sequence for ATAC).
    """
    # QC counters
    exact_match = nonexact_match = nonmatch = poly_g_barcode = poly_g_start = 0
    
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
            r1, q1, exact1 = check_putative_barcode(r1_str, r1_barcode_dict_exact, r1_barcode_dict_mismatch, q1_str)
            r2, q2, exact2 = check_putative_barcode(r2_str, r2_barcode_dict_exact, r2_barcode_dict_mismatch, q2_str)
            r3, q3, exact3 = check_putative_barcode(r3_str, r3_barcode_dict_exact, r3_barcode_dict_mismatch, q3_str)
            
            # check first ten base pairs of read 2 for homopolymer G
            if read_2.sequence[:10] == "G"*10:
                poly_g_start += 1
                
            # if corrected barcodes found, write to both R1 and R2 FASTQ files
            elif r1 and r2 and r3:
                if exact1 and exact2 and exact3:
                    exact_match += 1
                else:
                    nonexact_match += 1
                    
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
                    # create SequenceRecord object for read 2: use corrected header, 50bp read
                    corrected_read_2 = dnaio.SequenceRecord(corrected_header, read_2.sequence[:50], read_2.qualities[:50])
                
                # write to corrected FASTQ files
                writer.write(corrected_read_1, corrected_read_2)
                
            # check for homopolymer G in uncorrected barcode windows
            elif "G"*8 in r1_str and "G"*8 in r2_str and "G"*8 in r3_str:
                poly_g_barcode += 1
                
            else:
                nonmatch += 1
    
    # write QC stats
    with open(f"{prefix}_barcode_qc.txt", "w") as f:
        fields = ["library", "exact_match", "nonexact_match", "nonmatch", "poly_G_barcode", "poly_G_in_first_10bp"]
        f.write("\t".join(fields) + "\n")
        f.write("%s\t%s\t%s\t%s\t%s\t%s" % (prefix, exact_match, nonexact_match, nonmatch, poly_g_barcode, poly_g_start))
        

def main():
    args = parse_arguments()
    input_r1_fastq_file = getattr(args, "input_r1_fastq_file")
    input_r2_fastq_file = getattr(args, "input_r2_fastq_file")
    output_r1_fastq_file = getattr(args, "output_r1_fastq_file")
    output_r2_fastq_file = getattr(args, "output_r2_fastq_file")
    r1_barcode_file = getattr(args, "r1_barcode_file")
    r2_barcode_file = getattr(args, "r2_barcode_file")
    r3_barcode_file = getattr(args, "r3_barcode_file")
    sample_type = getattr(args, "sample_type")
    pkr = getattr(args, "pkr")
    prefix = getattr(args, "prefix")
    
    # read in barcodes, create dictionary including mismatches
    with open(r1_barcode_file) as f:
        (r1_barcode_dict_exact, r1_barcode_dict_mismatch) = create_barcode_dicts(f.read().strip().split())

    with open(r2_barcode_file) as f:
        (r2_barcode_dict_exact, r2_barcode_dict_mismatch) = create_barcode_dicts(f.read().strip().split())
        
    with open(r3_barcode_file) as f:
        (r3_barcode_dict_exact, r3_barcode_dict_mismatch) = create_barcode_dicts(f.read().strip().split())

    # write corrected FASTQs and QC stats
    process_fastqs(input_r1_fastq_file, input_r2_fastq_file,
                   output_r1_fastq_file, output_r2_fastq_file,
                   r1_barcode_dict_exact, r1_barcode_dict_mismatch, 
                   r2_barcode_dict_exact, r2_barcode_dict_mismatch, 
                   r3_barcode_dict_exact, r3_barcode_dict_mismatch,
                   sample_type, pkr, prefix)

if __name__ == "__main__":
    main()
