#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trim fastq
Removes dovetail (overlap) between R1 and R2
"""

import argparse
import dnaio
import Levenshtein

def parse_arguments():
    parser = argparse.ArgumentParser(description="Trim dovetail (overlap) between read1 and read2")
    parser.add_argument("input_r1_fastq_file", help="Filename for untrimmed input read 1 FASTQ file")
    parser.add_argument("input_r2_fastq_file", help="Filename for untrimmed input read 2 FASTQ file")
    parser.add_argument("output_r1_fastq_file", help="Filename for corrected output read 1 FASTQ file")
    parser.add_argument("output_r2_fastq_file", help="Filename for corrected output read 2 FASTQ file")
    
    return parser.parse_args()

# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}

def reverse_complement(sequence):
    """
    Return reverse complement of DNA sequence.
    """
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])

def process_fastqs(input_r1_fastq_file, input_r2_fastq_file,
                  output_r1_fastq_file, output_r2_fastq_file):

    # process FASTQs together
    with dnaio.open(input_r1_fastq_file, input_r2_fastq_file) as reader, \
         dnaio.open(output_r1_fastq_file, output_r2_fastq_file, mode="w") as writer:
             
        for read_1, read_2 in reader:
        # trim adapters for ATAC
            where = trim(read_1.sequence, read_2.sequence)
            if where > -1:
                # create SequenceRecord object for read 1 with trimmed sequence
                corrected_read_1 = dnaio.SequenceRecord(read_1.name, read_1.sequence[:where], read_1.qualities[:where])                
                # create SequenceRecord object for read 2 with trimmed sequence
                corrected_read_2 = dnaio.SequenceRecord(read_2.name, read_2.sequence[:where], read_2.qualities[:where])
            else:
                # no overlap, remains the same
                corrected_read_1 = read_1
                corrected_read_2 = read_2
            # write to corrected FASTQ files
            writer.write(corrected_read_1, corrected_read_2)

def trim(seq1,seq2):
    '''
    Find overlap between read1 and read2 and return location
    '''
    query = reverse_complement(seq2[0:20])
    idx = seq1.rfind(query) # look for perfect match
    if idx == -1:
        idx = fuzz_align(query,seq1)
    # found it, return everything through match
    if idx > -1:
        idx = idx+20
    else:
        idx = -1
    return idx

def fuzz_align(s_seq,l_seq):
    '''
    Align allowing Levenshtein distance of 1
    This iteration should go from the right end of l_seq
    since we want to do a rfind 
    '''
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq)
        if dist <= 1:  # find first then break
            return i
    return -1



def main():
    args = parse_arguments()
    input_r1_fastq_file = getattr(args, "input_r1_fastq_file")
    input_r2_fastq_file = getattr(args, "input_r2_fastq_file")
    output_r1_fastq_file = getattr(args, "output_r1_fastq_file")
    output_r2_fastq_file = getattr(args, "output_r2_fastq_file")
    process_fastqs(input_r1_fastq_file, input_r2_fastq_file,
                   output_r1_fastq_file, output_r2_fastq_file)


if __name__ == "__main__":
    main()
