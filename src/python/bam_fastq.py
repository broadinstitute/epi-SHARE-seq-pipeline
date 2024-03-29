#!/usr/bin/env python
"""
Write paired end reads from unmapped BAM file to FASTQ files.
Write only reads that match (given) barcodes within mismatch tolerance.
Output FASTQ files per sample based on splits of round 1 cell barcodes.

The BAM file contains RX and QX tags with sequence and quality, respectively, 
of the cell barcodes. The format is 10bp-10bp-9bp (dash separated) with the 
round 1 cell barcode in the first position, the round 2 cell barcode in the
second position, and the round 3 barcode in the third position. There is a 1 bp 
shift tolerance on either side of the barcode, except for round 3, where the
barcode is expected to be at the end. The expected location of the exact match 
is in the center.

The BAM file also contains an RG tag with the PKR identifier, so cells from the
same PKR (giving rise to separate ATAC and RNA libraries) may be matched.

The cell barcode and PKR id are added to the read name; for RNA, the UMI is 
added as well. The read2 FASTQ for RNA contains the cell barcode and UMI.

There are two FASTQ files written per sample as indicated by the round1 
barcode set file. Each line of the barcode set file represents a barcode set.
The first field is the name of the barcode set (e.g., 96-PLATE1-1st-THIRD) and 
the following fields are the barcodes. If samples are not split on the plate,
there will be only one line, with the first field as the name (e.g. 
96-FULL-PLATE) and the following fields the barcodes; in the case of one 96 
well plate, there will be 97 total fields on that one line. 

The round2 and round3 barcode files can be the same as the round1 barcode file
if the same barcodes are used for round2 and round3.

@author Neva Durand (c) 2021
"""

import argparse
import os
import sys
import pysam

# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}

def main(bam_file, r1_barcode_set_file, r2_barcode_file, r3_barcode_file, 
         sample_type, file_prefix):
    """
    Open involved files, read in barcode file and create dictionary,
    write unaligned BAM reads with barcodes matching 
    to FASTQ files.
    """
    # Read in barcodes, create dictionary including mismatches
    barcode_set, r1_barcodes = create_barcode_set(r1_barcode_set_file)
    (r1_barcode_dict_exact, 
        r1_barcode_dict_mismatch) = create_barcode_dict(r1_barcodes)
    with open(r2_barcode_file) as f:
        (r2_barcode_dict_exact, 
            r2_barcode_dict_mismatch) = create_barcode_dict(f.read().strip().split())
    with open(r3_barcode_file) as f:
        (r3_barcode_dict_exact, 
            r3_barcode_dict_mismatch) = create_barcode_dict(f.read().strip().split())

    barcode_set_values = set(barcode_set.values())

    # create file pointers from the round1 barcode set to write to
    left = dict()
    right = dict()
    for value in barcode_set_values:
        fp = open(file_prefix + '_' + value + '_R1.fastq', 'w')
        left[value] = fp
        fp = open(file_prefix + '_' + value + '_R2.fastq', 'w')
        right[value] = fp
    with pysam.Samfile(bam_file, 'rb', check_sq=False) as bam:
        process_bam(bam, left, right, r1_barcode_dict_exact, r1_barcode_dict_mismatch, 
                    r2_barcode_dict_exact, r2_barcode_dict_mismatch, 
                    r3_barcode_dict_exact, r3_barcode_dict_mismatch, 
                    barcode_set, sample_type, file_prefix)

def create_barcode_dict(barcode_list):
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

def create_barcode_set(file_path):
    """
    Create dictionary for looking up per round1 barcode which sample file
    to write to
    """    
    with open(file_path) as f:
        barcodeset = dict() # [seq: set_name]
        barcodelist = list() # [seq]
        for row in f.readlines():
            row = row.strip().split()
            name = row[0] # set name is first element
            for item in row[1:]:
                barcodeset[item] = name
                barcodelist.append(item)
    return barcodeset, barcodelist

def process_bam(bam, left, right, r1_barcode_dict_exact, r1_barcode_dict_mismatch, 
                r2_barcode_dict_exact, r2_barcode_dict_mismatch, 
                r3_barcode_dict_exact, r3_barcode_dict_mismatch, 
                barcode_set, sample_type, file_prefix):
    """
    Get reads from open BAM file and write them in pairs.
    left / right are arrays containing file pointers to left "read1" and right 
     "read2"
    The dictionaries are the round1, round2, round3 exact and single mismatch 
     hashtables
    The barcode set dictates which file (from left/right) is written to based 
     on the round1 barcode
    sample_type is RNA or ATAC
    file_prefix is for the QC stats output file; these stats indicate how many 
     reads are lost to various complexity issues (homopolymers in first 10bp of 
     read2, where UMI lives; homopolymers in cell barcode; cell barcode mismatch) 
    """
    
    qname = read_left = read_right = None
    exact = good = bad = poly_barcode = poly_umi = 0 # QC counters
    poly_g = 'GGGGGGGG'
    poly_g_umi = 'GGGGGGGGGG'

    for read in bam:
        if read.is_read1:
            # save and continue processing
            read_left = read
            qname = read_left.query_name
        else:
            # check that right/left query names are the same
            if qname == read.query_name:
                # check barcode via RX tag, change readname
                pkr_id = read.get_tag("RG")
                barcode_tag = read.get_tag("RX")
                R1str,R2str,R3str = barcode_tag.split("-",2)
                quality_tag = read.get_tag("QX")
                Q1str,Q2str,Q3str = quality_tag.split("-",2)
                R1 = R2 = R3 = None
                R1,Q1,exact1 = check_putative_barcode(R1str, r1_barcode_dict_exact, 
                                                      r1_barcode_dict_mismatch, Q1str)
                R2,Q2,exact2 = check_putative_barcode(R2str, r2_barcode_dict_exact, 
                                                      r2_barcode_dict_mismatch, Q2str)
                R3,Q3,exact3 = check_putative_barcode(R3str, r3_barcode_dict_exact, 
                                                      r3_barcode_dict_mismatch, Q3str)
                
                # examine "umi" region
                read_right = read
                umi = read_right.query_sequence[0:10]
                umi_qual = read_right.qual[0:10]

                if umi == poly_g_umi:
                    poly_umi += 1
                elif R1 and R2 and R3:
                    if exact1 and exact2 and exact3:
                        exact += 1
                    else:
                        good += 1
                        
                    # add cell barcodes to queryname
                    qname_barcode = ",".join([R1,R2,R3,pkr_id])
                    qname = qname + "_" + qname_barcode
                    
                    if sample_type == 'ATAC':
                        read_left.query_name = qname
                        read_right.query_name = qname
                        # left and right contain open file pointers
                        # write read to correct file based on R1 barcode
                        write_read(left[barcode_set[R1]], read_left)
                        write_read(right[barcode_set[R1]], read_right)
                    elif sample_type == 'RNA':
                        qname = qname + "_" + umi
                        read_left.query_name = qname
                        read_right.query_name = qname
                        read_right.query_sequence = R1 + R2 + R3 + umi
                        read_right.qual = Q1 + Q2 + Q3 + umi_qual
                        # left contains open file pointers
                        # write read to correct file based on R1 barcode
                        write_read(left[barcode_set[R1]], read_left)
                        write_read(right[barcode_set[R1]], read_right)
                elif poly_g in R1str and poly_g in R2str and poly_g in R3str:
                    poly_barcode += 1
                else:                 
                    bad += 1
    
    with open("qc.txt", 'w') as f:
        print("%s\t%s\t%s\t%s\t%s\t%s" % (file_prefix, exact, good, bad, poly_barcode, poly_umi), file=f)

def check_putative_barcode(barcode_str, barcode_dict_exact, barcode_dict_mismatch, quality_str):
    '''
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right
     shift
    '''
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


def write_read(fastq, read):
    """
    Write read to open FASTQ file.
    """
    info = {'index': int(not read.is_read1) + 1,
            'name':  read.query_name}
                        
    if read.is_reverse:
        info.update({'quality':  read.qual[::-1],
                     'sequence': reverse_complement(read.query_sequence)})
    else:
        info.update({'quality':  read.qual,
                     'sequence': read.query_sequence})
    fastq.write('@{name}\n{sequence}\n+\n{quality}\n'.format(**info))

def reverse_complement(sequence):
    """
    Return reverse complement of DNA sequence.
    """
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('bam_file', metavar='BAM_FILE',
                       help='file in BAM format to extract reads from')
    group.add_argument('r1_barcode_set_file',
                       help='file containing biosample splits in R1 barcodes, one split per line')
    group.add_argument('r2_barcode_file', 
                       help='file containing R2 barcodes, one line')
    group.add_argument('r3_barcode_file', 
                       help='file containing R3 barcodes, one line')
    group.add_argument('-p', dest='file_prefix', help='prefix for FASTQ files'
                       ' (default: BAM_FILE_R1.fq, BAM_FILE_R2.fq')
    group.add_argument('-s', dest='sample_type', help='sample type in this library'
                       ' (default: ATAC)', choices=['ATAC', 'RNA'], default='ATAC')
    args = parser.parse_args()

    # Default file_prefix
    name, _ = os.path.splitext(args.bam_file)
    file_prefix = args.file_prefix if args.file_prefix else name

    main(args.bam_file, args.r1_barcode_set_file, args.r2_barcode_file, 
         args.r3_barcode_file, args.sample_type, file_prefix)
