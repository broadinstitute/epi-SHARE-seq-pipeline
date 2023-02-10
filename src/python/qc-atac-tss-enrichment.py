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
def count_fragments_in_promoter(bam_filename,
                                tss_list,
                                flank=2000,
                                barcode_tag = "CB",
                                mapq_threshold = 30):
    """
    This funtion counts the per-barcode number of reads in the promoter region of the TSSs passed in input.
    The code follows the ArchR heuristic to minimize memory usage. Only the reads in the +/-50bp around the
    TSS will be counted. The control reagions are still the first and last 100bp of the interval.

    Parameters
    ----------
    bam_filename : str
        Path to the SAM/BAM file containing the mapped reads.
        File needs to be coordinate-sorted and indexed.
    tss_list : array
        Array containing the list of TSSs to be included.
        Each member of the array contains the following four elements:
        Chr, Start,  End, Strand
    flank : int
        How many base pairs to add upstream and downstream of the TSS.
        default: 2000
    barcode_tag : str
        Which tag in the BAM file contains the barcode id.
    mapq_threshold : int
        Keep only the reads with mapq score greater or equal.
        default: 30

    Returns
    -------
    Array
        Value: Aggregated number or reads per promoter. To produce the bulk TSS enrichment plot.

    Dictionary
        Key: Barcode
        Value: Number of reads around the TSS. (TSS+/-50)

    Dictionary
        Key: Barcode
        Value: Cumulative number of reads around the promoter (TSS +/-flank)
    """
    # TSS +/- flank
    promoter_size = flank * 2
    counts_dict_bulk = np.zeros(promoter_size)
    # To count the number of reads in the promoter region.
    fragments_in_promoter_counter = Counter()

    # 100(flank)signal TSS +/- 50(101)flank(100) = 301
    # Follwoing ArchR heuristic.
    counts_dict = defaultdict(lambda: np.zeros(301))

    bamfile = pysam.Samfile(bam_filename, "rb")

    for tss in tss_list:
        # TSS example: ["chr", "start", "end", "strand"]
        # Create the promoter region by adding the upstream and downstream.
        promoter_start = int(tss[1]) - flank
        promoter_end = int(tss[1]) + flank
        promoter_strand = tss[3]

        # Find all the fragments overlapping the promoter.
        for read in bamfile.fetch(str(tss[0]), max(0,promoter_start), promoter_end):
            # Check mapping quality
            if read.mapq < mapq_threshold or read.flag & 16 == 16:# or read.get_tag(barcode_tag)!="TCATCCTAGCAGAGCCTCGAGCGT":
                continue # Ignore low quality reads and reverse (coordinate-wise second) read in pair.

            fragment_start = read.pos + 4
            # The -1 is needed to get the last accessible element because BED/BAM files are 0-based [s,e)
            fragment_end = read.pos + read.tlen - 4
            fragment_length = fragment_end - fragment_start

            # Increment the counter for the specific barcode.
            barcode = read.get_tag(barcode_tag)
            fragments_in_promoter_counter[barcode] += 1

            # Update the array with the counts around the promoter.
            _add_read_to_dictionary(counts_dict_bulk,
                                    counts_dict,
                                    fragment_start,
                                    promoter_start,
                                    promoter_end,
                                    promoter_strand,
                                    barcode)

            _add_read_to_dictionary(counts_dict_bulk,
                                    counts_dict,
                                    fragment_end,
                                    promoter_start,
                                    promoter_end,
                                    promoter_strand,
                                    barcode)

    return counts_dict_bulk, counts_dict, fragments_in_promoter_counter

def _add_read_to_dictionary(fragment_counter,
                            fragment_counter_barcodes,
                            fragment_position,
                            promoter_start,
                            promoter_end,
                            promoter_strand,
                            barcode
                            ):
    weight = 1
    if  fragment_position >= promoter_start and fragment_position <= promoter_end - 1:
        # We are adding to an array that covers each basepair from promoter_start to promoter_end.
        # This converts the genomic coordinates to the index that we need to update in the array.
        index = int(fragment_position - promoter_start) # 0-based

        max_window = len(fragment_counter)-1 # max accessible array index
        tss_pos = int(max_window/2)

        add_count_to_barcode = False

        if index >= 0 and index <= 99:
            index_reduced = index
            add_count_to_barcode = True

        elif index >= max_window - 99:
            index_reduced = index - (max_window - 300)
            add_count_to_barcode = True

        elif index >= tss_pos - 50 and index <= tss_pos + 50:
            # Find the center of the region with max_window/2+50
            index_reduced = int(index - (tss_pos - 50)) + 100
            add_count_to_barcode = True

        # If the TSS is in the negative strand we need to add from the end.
        if promoter_strand == "-":
            index = max_window - index
            if add_count_to_barcode:
                index_reduced = 301 - index_reduced - 1

        # Finally update the dictionary
        fragment_counter[index] += weight
        if add_count_to_barcode:
            fragment_counter_barcodes[barcode][index_reduced] += weight
    return

def plot_tss_enrichment(raw_signal, smoothed_signal, out_file):
    fig=plt.figure(figsize=(8.0, 5.0))
    plt.plot(raw_signal,'k.')
    plt.plot(smoothed_signal,'r')
    #plt.xlabel('Position relative to center')
    plt.ylabel('TSS enrichment')
    plt.xticks([0, 2000, 4000], ['-2000', 'TSS', '+2000'])
    fig.savefig(out_file)
    plt.close(fig)

def compute_tss_enrichment(array_counts, window_size, png_file):
    size = len(array_counts)
    # We want to normalize the plot using the signal in the first and last 200 bps.
    normalization_factor = np.mean(array_counts[np.r_[0:100, size-100-1:size]])

    raw_signal = array_counts/max(0.2, normalization_factor)
    # Smooth using a window
    smoothed_signal = np.convolve(array_counts, np.ones(window_size),'same')/window_size/normalization_factor

    plot_tss_enrichment(raw_signal, smoothed_signal, png_file)

    return max(smoothed_signal)

def compute_tss_enrichment_barcode(array_counts):
    size = len(array_counts)
    # We want to normalize the plot using the signal in the first and last 200 bps.
    normalization_factor = 2*(np.sum(array_counts[np.r_[0:100, size-100:size]]))/200

    reads_in_tss = np.sum(array_counts[100:201])

    #normalization_factor = np.mean(array_counts[np.r_[0:100, 201:301]])
    #print(normalization_factor)
    # To avoid big TSS enrichment scores we use 0.2 as minimum. (Taken from ArchR).
    tss_enrichment = 2*reads_in_tss/101/max(0.2,normalization_factor)

    return reads_in_tss, tss_enrichment


if __name__ == '__main__':

    #args = _parse_sanitize_cmdline_arguments()

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("bam", help= "Path to the coordinate-sorted bam file.")
    parser.add_argument("-e", help= "Number of bases to extend to each side. (default= 1000)", type= int, default= 2000)
    parser.add_argument("-s", help="Column with strand information; 1-based. (default= 4)", type= int, default= 4)
    parser.add_argument("-w", "--window", help= "Smoothing window size for plotting. (default= 20)", type= int, default= 20)
    parser.add_argument("--bc_tag", help = "Specify the tag containing the cell barcode.", default="CB")
    parser.add_argument("--mapq_threshold", help= "Filter reads with a mapq value lower than the threshold.", type= int, default= 30)
    parser.add_argument("--prefix", help = "Prefix for the metrics output fil.")
    parser.add_argument("--tss", help= "TSS bed file")

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = args.bam[:-4]

    # Read arguments from command line
    args = parser.parse_args()

    # It is extremely fast. Don't think we need parallel processing.
    #cpus = len(os.sched_getaffinity(0))/2
    # Using column chr, start, end and what user input contains the strand information.
    tss_list = np.loadtxt(args.tss, 'str', usecols = (0,1,2,args.s-1))

    bulk_counts, barcode_counts, stats = count_fragments_in_promoter(args.bam,
                                                                     tss_list,
                                                                     flank= args.e,
                                                                     barcode_tag = args.bc_tag,
                                                                     mapq_threshold = args.mapq_threshold
                                                                    )

    per_barcode_output = f"{args.prefix}.tss_enrichment_barcode_stats.tsv"
    tss_enrichment_plot_fnp = f"{args.prefix}.tss_enrichment_bulk.png"

    with open(per_barcode_output,"w") as out_file:
        print(f"barcode\tfragments_in_promoter\treads_in_tss\ttss_enrichment", file=out_file)
        for barcode, fragments_in_promoter in stats.items():
            reads_in_tss, tss_enrichment = compute_tss_enrichment_barcode(barcode_counts[barcode])
            print(f"{barcode}\t{fragments_in_promoter}\t{reads_in_tss}\t{tss_enrichment}", file=out_file)

    with open(f"{args.prefix}.tss_score_bulk.txt", "w") as out_file:
        tss_score_bulk = compute_tss_enrichment(bulk_counts, args.window, tss_enrichment_plot_fnp)
        print(f"tss_Score", file=out_file)

