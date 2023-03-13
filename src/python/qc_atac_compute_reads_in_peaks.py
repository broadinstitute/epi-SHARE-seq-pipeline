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
def count_fragments_in_peaks(tabix_filename,
                                peaks_list,
                                mapq_threshold = 30):
    """
    This funtion counts the per-barcode number of reads in the peak region.

    Parameters
    ----------
    tabix_filename : str
        Path to the tabix file containing the fragments.
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
    reads_in_peaks_counter = defaultdict(set)
    fragments_in_peaks_counter = defaultdict(set)

    tabixfile = pysam.TabixFile(tabix_filename)

    for peak in peaks_list:
        peak_chr = str(peak[0])
        peak_start = int(peak[1])
        peak_end = int(peak[2])

        # Find all the fragments overlapping the promoter.
        for fragment in tabixfile.fetch(peak_chr, peak_start, peak_end):
            fragment_fields = fragment.split("\t")

            fragment_contig = fragment_fields[0]
            fragment_start = int(fragment_fields[1])
            fragment_end = int(fragment_fields[2])
            barcode = fragment_fields[3]

            fragment_id = "-".join(fragment_fields)
            fragments_in_peaks_counter[barcode].add(fragment_id)

            # Increment the counter for the specific barcode.
            if fragment_start >= peak_start and fragment_start <= peak_end-1:
                reads_in_peaks_counter[barcode].add(fragment_id+"start")

            if fragment_end >= peak_start and fragment_end <= peak_end-1:
                reads_in_peaks_counter[barcode].add(fragment_id+"end")

    return reads_in_peaks_counter, fragments_in_peaks_counter


if __name__ == '__main__':

    #args = _parse_sanitize_cmdline_arguments()

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("tabix", help= "Fragments file in tabix format and indexed.")
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

    reads_in_peaks, fragments_in_peaks = count_fragments_in_peaks(args.tabix,
                                                  peaks_list
                                                 )
    output_fnp = f"{prefix}.reads.in.peak.tsv"

    with open(output_fnp,"w") as out_file:
        print(f"barcode\treads_peaks\tfragment_peaks", file=out_file)
        for barcode,fragments_in_peak in fragments_in_peaks.items():
            print(f"{barcode}\t{len(reads_in_peaks[barcode])}\t{len(fragments_in_peak)}", file=out_file)
