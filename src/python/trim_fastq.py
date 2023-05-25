#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trim fastq
Removes dovetail (overlap) between R1 and R2
"""

import argparse
import Levenshtein
import xopen
from collections import deque

def parse_arguments():
    parser = argparse.ArgumentParser(description="Trim dovetail (overlap) between read1 and read2")
    parser.add_argument("input_read1_fastq_file", help="Filename for untrimmed input read 1 FASTQ file")
    parser.add_argument("input_read2_fastq_file", help="Filename for untrimmed input read 2 FASTQ file")
    parser.add_argument("output_read1_fastq_file", help="Filename for corrected output read 1 FASTQ file")
    parser.add_argument("output_read2_fastq_file", help="Filename for corrected output read 2 FASTQ file")
    parser.add_argument("trimming_stats_file", help="Filename for txt file containing trimming statistics")
    
    return parser.parse_args()

REV_COMP = str.maketrans("ATGC", "TACG")
def reverse_complement(seq):
    return str.translate(seq, REV_COMP)[::-1]

def trim_fastqs(input_read1_fastq_file, input_read2_fastq_file,
                output_read1_fastq_file, output_read2_fastq_file,
                trimming_stats_file):
    """
    Trim reads if overlapping, write reads to output FASTQ files.
    Produces file enumerating how many reads were processed and trimmed.
    """
    # counters
    total = trimmed = 0
    
    read1_out_writer = xopen.xopen(output_read1_fastq_file, mode="w")
    read2_out_writer = xopen.xopen(output_read2_fastq_file, mode="w")

    buffer1 = deque()
    buffer2 = deque()
    buffer_counter = 0
    
    # process FASTQs together
    with xopen.xopen(input_read1_fastq_file, mode= "r", threads= 8) as read1_fh, xopen.xopen(input_read2_fastq_file, mode= "r", threads= 8) as read2_fh:
        for readline1, readline2 in zip(read1_fh, read2_fh):
            total += 2
            
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
             
            # trim adapters for ATAC
            where = trim(sequence1, sequence2)
            
            if where > -1:
                trimmed += 2
                
                # add trimmed read 1 to buffer
                trimmed_read1 = f"{name1}\n{sequence1[:where]}\n+\n{quality1[:where]}\n"
                buffer1.append(trimmed_read1)
                
                # add trimmed read 2 to buffer
                trimmed_read2 = f"{name2}\n{sequence2[:where]}\n+\n{quality2[:where]}\n"
                buffer2.append(trimmed_read2)
                
            else:
                # add original read 1 to buffer
                read1 = f"{name1}\n{sequence1}\n+\n{quality1}\n"
                buffer1.append(read1)
                
                # add original read 1 to buffer
                read2 = f"{name2}\n{sequence2}\n+\n{quality2}\n"
                buffer2.append(read2)
                
            buffer_counter += 1    
                
            # write reads to trimmed FASTQ files 
            if buffer_counter == 10000000:
                read1_out_writer.write("".join(buffer1))
                buffer1.clear()
                read2_out_writer.write("".join(buffer2))
                buffer2.clear()
                buffer_counter = 0
    
    # write out remaining reads
    if buffer_counter > 0:
        read1_out_writer.write("".join(buffer1))
        buffer1.clear()
        read2_out_writer.write("".join(buffer2))
        buffer2.clear()
        buffer_counter = 0
    
    # write trimming statistics output file
    with open(trimming_stats_file, "w") as f:
        fields = ["total_reads", "untrimmed_reads", "trimmed_reads", "%trimmed"]
        f.write("\t".join(fields) + "\n")
        f.write("%i\t%i\t%i\t%0.1f" % (total, total-trimmed, trimmed, trimmed/total*100))

def trim(seq1, seq2):
    """
    Find overlap between read1 and read2 and return location
    """
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

def fuzz_align(s_seq, l_seq):
    """
    Align allowing Levenshtein distance of 1
    This iteration should go from the right end of l_seq
    since we want to do a rfind 
    """
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq, score_cutoff= 1)
        if dist <= 1:  # find first then break
            return i
    return -1

def main():
    args = parse_arguments()
    input_read1_fastq_file = getattr(args, "input_read1_fastq_file")
    input_read2_fastq_file = getattr(args, "input_read2_fastq_file")
    output_read1_fastq_file = getattr(args, "output_read1_fastq_file")
    output_read2_fastq_file = getattr(args, "output_read2_fastq_file")
    trimming_stats_file = getattr(args, "trimming_stats_file")
    
    trim_fastqs(input_read1_fastq_file, input_read2_fastq_file,
                output_read1_fastq_file, output_read2_fastq_file,
                trimming_stats_file)


if __name__ == "__main__":
    main()
