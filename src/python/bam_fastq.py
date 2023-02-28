#!/usr/bin/env python
"""
Write paired end reads from unmapped BAM file to FASTQ files.
Write only reads that match (given) barcodes within 
mismatch tolerance.

@author Neva Durand (c) 2021
"""

import argparse
import Levenshtein
import os
import sys
import pysam

# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}

def main(bam_file, r1_barcode_sets, r2_barcode_file, r3_barcode_file, pkr_id, sample_type, file_prefix=None):
    """
    Open involved files, read in barcode file and create dictionary,
    write unaligned BAM reads with barcodes matching 
    to FASTQ files.
    """
    barcode_set, r1_barcodes = create_barcode_set(r1_barcode_sets)
    r1_barcode_dict = create_barcode_dict(r1_barcodes)
    # Read in barcodes, create dictionary including mismatches
    with open(r2_barcode_file) as f:
        r2_barcode_dict = create_barcode_dict(f.read().strip().split())
    with open(r3_barcode_file) as f:
        r3_barcode_dict = create_barcode_dict(f.read().strip().split())

    barcode_set_values = set(barcode_set.values())
    name, _ = os.path.splitext(bam_file)
    if not file_prefix:
        file_prefix = name
    fnames = []
    left = dict()
    right = dict()
    for value in barcode_set_values:
        fp = open(file_prefix + '_' + value + '_R1.fastq', 'w')
        left[value] = fp
        fp = open(file_prefix + '_' + value + '_R2.fastq', 'w')
        right[value] = fp
    with pysam.Samfile(bam_file, 'rb', check_sq=False) as bam:
        process_bam(bam, left, right, r1_barcode_dict, r2_barcode_dict, r3_barcode_dict, barcode_set, pkr_id, sample_type, file_prefix)

def create_barcode_dict(barcode_list):
    """
    Adds each barcode and its mismatch possibilities to the dictionary
    """
    barcode_dict = dict()
    for code in barcode_list:
        barcode_dict[code]=code # exact match
        for i, c in enumerate(code):
            for base in 'ACGTN':
                if c != base:
                    # add mismatch possibilities at pos i
                    barcode_dict[(code[:i] + base + code[i + 1:])] = code
    return barcode_dict    

def create_barcode_set(file_path):
    with open(file_path) as f:
        barcodeset = dict()
        barcodelist = list()
        for row in f.readlines():
            row = row.strip().split()
            name = row[0]
            for item in row[1:]:
                barcodeset[item] = name
                barcodelist.append(item)
    return barcodeset, barcodelist

def process_bam(bam, left, right, r1_barcode_dict, r2_barcode_dict, r3_barcode_dict, barcode_set, pkr_id, sample_type, file_prefix):
    """
    Get reads from open BAM file and write them in pairs.
    """
    
    qname = read_left = read_right = None
    good = bad = ggg = homop = 0
    G = 'GGGGGGGG'

    for read in bam:
        if read.is_read1:
            # save and continue processing
            read_left = read
            qname = read_left.query_name
        else:
            # check that right/left query names are the same
            if qname == read.query_name:
                # check barcode via RX tag, change readname
                barcode_tag = read.get_tag("RX")
                R1str,R2str,R3str = barcode_tag.split("-",2)
                quality_tag = read.get_tag("QX")
                Q1str,Q2str,Q3str = quality_tag.split("-",2)
                R1 = R2 = R3 = None
                R1,Q1 = check_putative_barcode(R1str, r1_barcode_dict, Q1str)
                R2,Q2 = check_putative_barcode(R2str, r2_barcode_dict, Q2str)
                R3,Q3 = check_putative_barcode(R3str, r3_barcode_dict, Q3str)
                
                if R1 and R2 and R3:
                    good += 1
                    # add cell barcodes to queryname
                    qname_barcode = ",".join([R1,R2,R3,pkr_id])
                    qname = qname + "_" + qname_barcode
                    read_right = read
                    if sample_type == 'ATAC':
                        read_left.query_name = qname
                        read_right.query_name = qname
                        # trim adapters for ATAC
                        where = trim(read_left.query_sequence, read_right.query_sequence)
                        # left and right contain open file pointers
                        # write read to correct file based on R1 barcode
                        write_read(left[barcode_set[R1]], read_left, where)
                        write_read(right[barcode_set[R1]], read_right, where)
                    elif sample_type == 'RNA':
                        # add UMI to queryname
                        umi = read_right.query_sequence[0:10]
                        umi_qual = read_right.qual[0:10]
                        #umi_qual = ''.join(map(lambda x: chr( x+33 ), read_right.query_qualities[0:10]))
                        qname = qname + "_" + umi
                        read_left.query_name = qname
                        read_right.query_name = qname
                        read_right.query_sequence = R1 + R2 + R3 + umi
                        read_right.qual = Q1 + Q2 + Q3 + umi_qual
                        # left contains open file pointers
                        # write read to correct file based on R1 barcode
                        write_read(left[barcode_set[R1]], read_left)
                        write_read(right[barcode_set[R1]], read_right)
                elif G in R1str and G in R2str and G in R3str:
                    read_right = read
                    umi = read_right.query_sequence[0:10]
                    if umi == 'GGGGGGGGGG':
                        homop += 1
                    else:
                        ggg += 1
                else:                 
                    bad += 1
    
    with open("qc.txt", 'w') as f:
        print("%s\t%s\t%s\t%s\t%s" % (file_prefix, good, bad, ggg, homop), file=f)

def check_putative_barcode(barcode_str, barcode_dict, quality_str):
    '''
    Procedure: check exact match R1, the 1bp left/right shift
    In the future the exact match will allow 2 mismatches
    and the 1bp left/right only 1 mismatch
    '''
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
    return value,quality


def trim(seq1,seq2):
    '''
    Trim putative adapters in ATAC reads
    This code is buggy (idx > 0) and ought to be verified
    Not clear it's better than more standard adapter trimming
    Need to check alignment stats
    Also Levenshtein vs Hamming, speed & accuracy
    '''
    query = reverse_complement(seq2[0:20])
    idx = seq1.rfind(query) # look for perfect match
    if idx == -1:
        idx = fuzz_align(query,seq1)
    # found it, return everything through match
    # NOTE: idx > 0 is incorrect
    if idx > 0:
        idx = idx+20-1
    else:
        idx = -1
    return idx

def fuzz_align(s_seq,l_seq):
    '''
    Check tradeoff using Hamming instead of Levenshtein
    This iteration should go from the right end of l_seq
    since we want to do a rfind 
    '''
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq)
        if dist <= 1:  # find first then break
            return i
    return -1

def write_read(fastq, read, idx=-1):
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
    if idx > -1:
        info.update({'quality':  read.qual[0:idx],
                    'sequence': read.query_sequence[0:idx]})
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
    group.add_argument('r1_barcode_sets',
                       help='file containing biosample splits in R1 barcodes, one split per line')
    group.add_argument('r2_barcode_file', 
                       help='file containing R2 barcodes, one line')
    group.add_argument('r3_barcode_file', 
                       help='file containing R3 barcodes, one line')
    group.add_argument('pkr_id',
                       help='PKR ID, should be matched ATAC and RNA')
    group.add_argument('-p', dest='file_prefix', help='prefix for FASTQ files'
                       ' (default: BAM_FILE_R1.fq, BAM_FILE_R2.fq')
    group.add_argument('-s', dest='sample_type', help='sample type in this library'
                       ' (default: ATAC)', choices=['ATAC', 'RNA'], default='ATAC')
    args = parser.parse_args()
    main(args.bam_file, args.r1_barcode_sets, args.r2_barcode_file, args.r3_barcode_file, args.pkr_id, args.sample_type, args.file_prefix)
