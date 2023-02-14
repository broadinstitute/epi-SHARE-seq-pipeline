#!/usr/bin/env python3

# Author: Eugenio Mattei, Broad Institute of MIT and Harvard
# modified from Jason Buenrostro's tool

import argparse
import os
import pysam
from collections import Counter
from collections import defaultdict
import numpy as np

#import os
#import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from multiprocessing import Pool


##### DEFINE FUNCTIONS #####
def count_fragments_in_peaks(bam_filename,
                                peaks_list,
                                barcode_tag = "CB",
                                mapq_threshold = 30):
    """
    This funtion counts the per-barcode number of reads in the peak region.

    Parameters
    ----------
    bam_filename : str
        Path to the SAM/BAM file containing the mapped reads.
        File needs to be coordinate-sorted and indexed.
    peaks_list : array
        Array containing the list of peaks to be included.
        Each member of the array contains the following four elements:
        Chr, Start,  End, Strand
    barcode_tag : str
        Which tag in the BAM file contains the barcode id.
    mapq_threshold : int
        Keep only the reads with mapq score greater or equal.
        default: 30

    Returns
    -------

    Dictionary
        Key: Barcode
        Value: Number of fragments in peaks.
    """
    # To count the number of fragments in peaks
    fragments_in_peaks_counter = defaultdict(set)

    bamfile = pysam.Samfile(bam_filename, "rb")

    for peak in peaks_list:
        peak_chr = str(peak[0])
        peak_start = int(peak[1])
        peak_end = int(peak[2])

        # Find all the fragments overlapping the promoter.
        for read in bamfile.fetch(peak_chr, peak_start, peak_end):
            # Check mapping quality
            if read.mapq < mapq_threshold or read.flag & 16 == 16:# or read.get_tag(barcode_tag)!="TCATCCTAGCAGAGCCTCGAGCGT":
                continue # Ignore low quality reads and reverse (coordinate-wise second) read in pair.

            # Increment the counter for the specific barcode.
            barcode = read.get_tag(barcode_tag)
            fragments_in_peaks_counter[barcode].add(read.name)

    return fragments_in_peaks_counter


if __name__ == '__main__':

    #args = _parse_sanitize_cmdline_arguments()

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("bam", help= "Path to the coordinate-sorted bam file.")
    parser.add_argument("--bc_tag", help = "Specify the tag containing the cell barcode.", default="CB")
    parser.add_argument("--mapq_threshold", help= "Filter reads with a mapq value lower than the threshold.", type= int, default= 30)
    parser.add_argument("--prefix", help = "Prefix for the metrics output fil.")
    parser.add_argument("--peaks", help= "Peaks bed file")

    # Read arguments from command line
    args = parser.parse_args()

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = args.bam[:-4]

    # It is extremely fast. Don't think we need parallel processing.
    #cpus = len(os.sched_getaffinity(0))/2
    # Using column chr, start, end and what user input contains the strand information.
    peaks_list = np.loadtxt(args.peaks, 'str', usecols = (0,1,2))

    fragments_in_peaks = count_fragments_in_peaks(args.bam,
                                                     peaks_list,
                                                     barcode_tag = args.bc_tag,
                                                     mapq_threshold = args.mapq_threshold
                                                    )
    output_fnp = f"{prefix}.fragments.in.peak.tsv"

    with open(output_fnp,"w") as out_file:
        print(f"barcode\tfragments_peaks", file=out_file)
        for barcode,fragments_in_peak in fragments_in_peaks.items():
            print(f"{barcode}\t{len(fragments_in_peak)}", file=out_file)
