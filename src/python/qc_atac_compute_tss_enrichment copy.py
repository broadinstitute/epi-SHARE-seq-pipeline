#!/usr/bin/env python3

# Author: Eugenio Mattei, Broad Institute of MIT and Harvard
# modified from Jason Buenrostro's tool

import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path as path
import pysam
import xopen
from collections import defaultdict
from tempfile import mkdtemp
#from multiprocessing import Pool

matplotlib.use('Agg')

##### DEFINE FUNCTIONS #####
def count_fragments_in_promoter(tabix_filename,
                                tss_list,
                                flank=2000):
    """
    This funtion counts the per-barcode number of reads in the promoter region of the TSSs passed in input.
    The code follows the ArchR heuristic to minimize memory usage. Only the reads in the +/-50bp around the
    TSS will be counted. The control reagions are still the first and last 100bp of the interval.

    Parameters
    ----------
    tabix_filename : str
        Path to the tabix file containing the fragments.
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

    # Get total number of barcodes
    unique_barcodes = set()
    # this will be used to access the correct row in the numpy matrix.
    barcodes_index_map = dict()
    total_barcodes = 0

    # Unique the barcodes
    tabix_fh = xopen.xopen(tabix_filename, mode="r", threads=4)
    for fragment in pysam.tabix_iterator(tabix_fh, pysam.asTuple()):
        unique_barcodes.add(fragment[3])
        
    for bc in unique_barcodes:
        barcodes_index_map[bc] = total_barcodes
        total_barcodes += 1
    
    # Cleaning up
    tabix_fh.close()
    del unique_barcodes

    # TSS +/- flank
    promoter_size = flank * 2
    counts_dict_bulk = np.zeros(promoter_size)

    # To count the number of reads in the promoter region.
    fragments_in_promoter_counter = defaultdict(set)
    reads_in_tss_counter = defaultdict(set)
    reads_in_promoter_counter = defaultdict(set)

    # 100(flank)signal TSS +/- 50(101)flank(100) = 301
    # Follwoing ArchR heuristic.
    counts_dict = defaultdict(lambda: np.zeros(301))
    print("created file for memmap - file.data")
    count_matrix = np.memmap('file.data', mode="w+",
                             dtype=np.int16,
                             shape=(total_barcodes, 301))

    tabixfile = pysam.TabixFile(tabix_filename)
    for tss in tss_list:
        # print(tss)
        # TSS example: ["chr", "start", "end", "strand"]
        # Create the promoter region by adding the upstream and downstream.
        promoter_start = int(tss[1]) - flank
        promoter_end = int(tss[1]) + flank
        promoter_strand = tss[3]
        counter = 0
        # if not(tss[0] == "chr1" and tss[1] == "24614812" and tss[2] == "24614813" and tss[3] == "-"):
        #     continue
        # print(f"{promoter_start}\t{promoter_end}\t{promoter_strand}")
        # Find all the fragments overlapping the promoter.
        try:
            for fragment in tabixfile.fetch(str(tss[0]),
                                            max(0, promoter_start),
                                            promoter_end):

                fragment_fields = fragment.split("\t")

                # fragment_contig = fragment_fields[0]
                fragment_start = int(fragment_fields[1])
                fragment_end = int(fragment_fields[2])
                barcode = fragment_fields[3]

                fragment_id = "-".join(fragment_fields)

                counter += 1

                # fragment_length = fragment_end - fragment_start

                # Increment the counter for the specific barcode.
                fragments_in_promoter_counter[barcode].add(fragment_id)

                # Update the array with the counts around the promoter.
                # I need 'start' and 'end' to create a unique id for the reads and I can count
                # the reads in TSSs wihtout double counting.
                _add_read_to_dictionary(counts_dict_bulk,
                                        count_matrix,
                                        fragment_start,
                                        promoter_start,
                                        promoter_end,
                                        promoter_strand,
                                        barcode,
                                        barcodes_index_map[barcode],
                                        reads_in_tss_counter,
                                        reads_in_promoter_counter,
                                        fragment_id,
                                        "start")

                _add_read_to_dictionary(counts_dict_bulk,
                                        count_matrix,
                                        fragment_end,
                                        promoter_start,
                                        promoter_end,
                                        promoter_strand,
                                        barcode,
                                        barcodes_index_map[barcode],
                                        reads_in_tss_counter,
                                        reads_in_promoter_counter,
                                        fragment_id,
                                        "end")
                #count_matrix.flush()
        except ValueError:
            print(f"No reads found for {tss[0]}:{max(1,promoter_start)}-{promoter_end}.")

    return barcodes_index_map, counts_dict_bulk, count_matrix, fragments_in_promoter_counter, reads_in_tss_counter, reads_in_promoter_counter

def _add_read_to_dictionary(fragment_counter,
                            fragment_counter_barcodes,
                            fragment_position,
                            promoter_start,
                            promoter_end,
                            promoter_strand,
                            barcode,
                            barcode_index,
                            reads_in_tss,
                            reads_in_promoter,
                            fragment_id,
                            end
                            ):
    
    weight = 1
    if fragment_position >= promoter_start and \
       fragment_position <= promoter_end - 1:
        
        # We are adding to an array that covers each basepair from promoter_start to promoter_end.
        # This converts the genomic coordinates to the index that we need to update in the array.
        index = int(fragment_position - promoter_start) # 0-based
        index_reduced = -1

        max_window = len(fragment_counter) - 1  # max accessible array index
        tss_pos = int(max_window/2)

        add_count_to_barcode = False
        add_count_to_barcode_tss = False

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
            add_count_to_barcode_tss = True

        # If the TSS is in the negative strand we need to add from the end.
        if promoter_strand == "-":
            index = max_window - index
            if add_count_to_barcode:
                index_reduced = 301 - index_reduced - 1

        # Finally update the dictionary
        fragment_counter[index] += weight  # bulk tss enrichment
        reads_in_promoter[barcode].add(fragment_id+end)

        if add_count_to_barcode:
            fragment_counter_barcodes[barcode_index][index_reduced] += weight  # per barcode tss enrichment
        if add_count_to_barcode_tss:
            reads_in_tss[barcode].add(fragment_id)
        
    return


def plot_tss_enrichment(raw_signal, smoothed_signal, out_file):
    
    fig = plt.figure(figsize=(8.0, 5.0))
    plt.plot(raw_signal, 'k.')
    plt.plot(smoothed_signal, 'r')
#   plt.xlabel('Position relative to center')
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
    smoothed_signal = np.convolve(array_counts, np.ones(window_size), 'same')/window_size/normalization_factor

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
    tss_enrichment = 2*reads_in_tss/101/max(0.2, normalization_factor)

    return tss_enrichment, reads_in_tss


if __name__ == '__main__':

    #args = _parse_sanitize_cmdline_arguments()

    msg = "Add the description"
    parser = argparse.ArgumentParser(description=msg)

    # Adding optional argument
    parser.add_argument("tabix", help="Fragments file in tabix format and indexed.")
    parser.add_argument("-e", help="Number of bases to extend to each side. (default= 2000)", type=int, default=2000)
    parser.add_argument("-s", help="Column with strand information; 1-based. (default= 4)", type=int, default=4)
    parser.add_argument("-w", "--window", help="Smoothing window size for plotting. (default= 20)", type=int, default=20)
    parser.add_argument("--prefix", help="Prefix for the metrics output file.")
    parser.add_argument("--tss", help="TSS bed file")

    # Read arguments from command line
    args = parser.parse_args()

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = args.bam[:-4]

    # It is extremely fast. Don't think we need parallel processing.
    #cpus = len(os.sched_getaffinity(0))/2
    # Using column chr, start, end and what user input contains the strand information.
    tss_list = np.loadtxt(args.tss, 'str', usecols=(0, 1, 2, args.s-1))

    barcode_index_map, bulk_counts, barcode_counts, stats, reads_in_tss, reads_in_promoter = count_fragments_in_promoter(args.tabix,
                                                                                                      tss_list,
                                                                                                      flank=args.e,
                                                                                                      )

    per_barcode_output = f"{args.prefix}.tss_enrichment_barcode_stats.tsv"
    tss_enrichment_plot_fnp = f"{args.prefix}.tss_enrichment_bulk.png"

    with open(per_barcode_output, "w") as out_file:
        # reads_tss_total exists to match the counts produced by ArchR.
        # A read can be counted twice because in two different TSSs.
        # I report the number of reads in TSS without double counting
        # but also the total count like ArchR
        print(f"barcode\tfragments_promoter\treads_tss\treads_promoter\ttss_enrichment\treads_tss_total", file=out_file)
        for barcode, fragments_in_promoter in stats.items():
            tss_enrichment, reads_sum = compute_tss_enrichment_barcode(barcode_counts[barcode_index_map[barcode], :])
            print(f"{barcode}\t{len(fragments_in_promoter)}\t{len(reads_in_tss[barcode])}\t{len(reads_in_promoter[barcode])}\t{tss_enrichment}\t{reads_sum}", file=out_file)

    with open(f"{args.prefix}.tss_score_bulk.txt", "w") as out_file:
        tss_score_bulk = compute_tss_enrichment(bulk_counts, args.window, tss_enrichment_plot_fnp)
        print(f"tss_enrichment\n{tss_score_bulk}", file=out_file)
