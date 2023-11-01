#!/usr/bin/env python3

# Author: Eugenio Mattei, Broad Institute of MIT and Harvard
# modified from Jason Buenrostro's tool

import argparse
import itertools
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import operator
import pysam

from collections import defaultdict
from functools import partial
from typing import BinaryIO, TextIO

matplotlib.use('Agg')


def get_outputs(idx, bulk_signal_counter, barcode_signal_counter, barcode_statistics):
    return {
        "idx": idx,
        "bulk_signal": bulk_signal_counter,
        "barcode_signal":  barcode_signal_counter,
        "barcode_statistics": barcode_statistics
        }

def compute_tss_score(idx: int,
                      tabix_filename: BinaryIO,
                      regions_bed: TextIO,
                      flank: int = 2000):
    """
    This funtion counts the per-barcode number of reads in the promoter region of the TSSs passed in input.
    The code follows the ArchR heuristic to minimize memory usage. Only the reads in the +/-50bp around the
    TSS will be counted. The control reagions are still the first and last 100bp of the interval.

    Parameters
    ----------
    tabix_filename : str
        Path to the tabix file containing the fragments.
        File needs to be coordinate-sorted and indexed.
    regions_bed : array
        Array containing the list of regions to be included.
        Each member of the array contains the following four elements:
        Chr, Start,  End, Strand
    flank : int
        How many base pairs to add upstream and downstream of the TSS.
        default: 2000

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
    bulk_signal_counter = np.zeros(promoter_size)
    # 100(flank)signal TSS +/- 50(101)flank(100) = 301
    # Follwoing ArchR heuristic.
    barcode_signal_counter = defaultdict(partial(np.zeros, 301))
    # To count the number of reads in the promoter region.
    # Array: [0]: reads in tss
    #        [1]: reads in promoter
    #        [2]: fragments in promoter
    barcode_statistics = defaultdict(partial(np.zeros,3))
    tabixfile = pysam.TabixFile(tabix_filename)
    for tss in regions_bed:
        # Create the promoter region by adding the upstream and downstream.
        promoter_start = int(tss[1]) - flank
        promoter_end = int(tss[1]) + flank
        promoter_strand = tss[3]
        try:
            for fragment in tabixfile.fetch(str(tss[0]), max(0, promoter_start), promoter_end):
                fragment_fields = fragment.split("\t")

#                fragment_contig = fragment_fields[0]
                fragment_start = int(fragment_fields[1])
                fragment_end = int(fragment_fields[2])
                barcode = fragment_fields[3]

                # I am using the helper so I don't have to
                # access the dictionary constantly.
                barcode_stats_helper = barcode_statistics[barcode]

                # Even if the fragment overlaps multiple promoter
                # I would like to count it once. Computationally
                # expensive so I gave up. ArchR does the same.
#                fragment_id = "-".join(fragment_fields)
#                fragment_length = fragment_end - fragment_start

                # Fragments in promoter counter increased for
                # the barcode.
                barcode_stats_helper[2] += 1

                # Update the array with the counts around the promoter.
                # I need 'start' and 'end' to create a unique id for the reads and I can count
                # the reads in TSSs wihtout double counting.
                _add_read_to_dictionary(bulk_signal_counter,
                                        barcode_signal_counter,
                                        fragment_start,
                                        promoter_start,
                                        promoter_end,
                                        promoter_strand,
                                        barcode,
                                        barcode_stats_helper)

                _add_read_to_dictionary(bulk_signal_counter,
                                        barcode_signal_counter,
                                        fragment_end,
                                        promoter_start,
                                        promoter_end,
                                        promoter_strand,
                                        barcode,
                                        barcode_stats_helper)
                # Saving the statistics per barcode back.
                barcode_statistics[barcode] = barcode_stats_helper
        except ValueError:
            print(f"No reads found for {tss[0]}:{max(1,promoter_start)}-{promoter_end}.")

    return get_outputs(idx, bulk_signal_counter, barcode_signal_counter, barcode_statistics)


def _add_read_to_dictionary(bulk_signal_counter,
                            barcode_signal_counter,
                            fragment_position,
                            promoter_start,
                            promoter_end,
                            promoter_strand,
                            barcode,
                            stats_counter):
    if fragment_position < promoter_start or fragment_position > promoter_end - 1:
        # This position does not cover the promoter.
        return
    # We are adding to an array that covers each basepair from promoter_start to promoter_end.
    # This converts the genomic coordinates to the index that we need to update in the array.
    index = int(fragment_position - promoter_start)  # 0-based
    index_reduced = -1

    max_window = len(bulk_signal_counter) - 1  # max accessible array index
    tss_pos = int(max_window/2)
    weight = 0

    if index >= 0 and index <= 99:
        index_reduced = index
        weight = 1

    elif index >= max_window - 99:
        index_reduced = index - (max_window - 300)
        weight = 1

    elif index >= tss_pos - 50 and index <= tss_pos + 50:
        # Find the center of the region with max_window/2+50
        index_reduced = int(index - (tss_pos - 50)) + 100
        weight = 1
        stats_counter[0] += 1

    # If the TSS is in the negative strand we need to add from the end.
    if promoter_strand == "-":
        index = max_window - index
        if weight > 0:
            index_reduced = 301 - index_reduced - 1

    # Finally update the dictionary
    bulk_signal_counter[index] += 1  # Bulk TSS enrichment
    stats_counter[1] += 1  # Reads in TSS

    barcode_signal_counter[barcode][index_reduced] += weight  # Per barcode TSS enrichment
    return


def plot_tss_enrichment(raw_signal, smoothed_signal, out_file):
    fig = plt.figure(figsize=(8.0, 5.0))
    plt.plot(raw_signal, 'k.')
    plt.plot(smoothed_signal, 'r')
#    plt.xlabel('Position relative to center')
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
#    normalization_factor = np.mean(array_counts[np.r_[0:100, 201:301]])
    # To avoid big TSS enrichment scores we use 0.2 as minimum. (Taken from ArchR).
    tss_enrichment = 2*reads_in_tss/101/max(0.2, normalization_factor)

    return tss_enrichment, reads_in_tss


def _prepare_args_for_counting(
        tabix_file: BinaryIO,
        regions_list: list,
        flank_size: int,
        cpus: int):
    # Array storing the list of arguments
    # for each thread.
    input_args = []
    # Create a number of chunks equal to
    # the number of cpus passed in input.
    for idx, chunk in enumerate(np.array_split(regions_list, cpus)):
        input_args.append(
            (idx,
             tabix_file,
             chunk,
             flank_size)
            )

    return input_args


def _merge_barcode_results(barcode, results, barcode_stats):
    tmp_signal = np.zeros(301)
    tmp_statistic = np.zeros(3)

    for result, statistic in zip(results, barcode_stats):
        if barcode in result:
            tmp_signal += result[barcode]
        if barcode in statistic:
            tmp_statistic += statistic[barcode]

    return tmp_signal, tmp_statistic


def _merge_bulk_signal(results):
    tmp = np.ones(len(results[0]))
    for result in results:
        tmp += result
    return tmp


if __name__ == '__main__':

    # args = _parse_sanitize_cmdline_arguments()

    msg = "Add the description"
    parser = argparse.ArgumentParser(description=msg)

    # Adding optional argument
    parser.add_argument("tabix", help="Fragments file in tabix format and indexed.")
    parser.add_argument("-e", help="Number of bases to extend to each side. (default= 2000)", type=int, default=2000)
    parser.add_argument("-s", help="Column with strand information; 1-based. (default= 4)", type=int, default=4)
    parser.add_argument("-p", help="Number or threads for parallel processing (default= 1)", type=int, default=1)
    parser.add_argument("-w", "--window",
                        type=int, default=20, help="Smoothing window size for plotting. (default= 20)")
    parser.add_argument("--prefix", help="Prefix for the metrics output file.")
    parser.add_argument("--regions", help="Bed file with the regions of interest")

    # Read arguments from command line
    args = parser.parse_args()

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = "sample"

    # One possible way to get the number of cores.
#    cpus = len(os.sched_getaffinity(0))/2
    # Using column chr, start, end and what user input contains the strand information.
    regions_list = np.loadtxt(args.regions, 'str', usecols=(0, 1, 2, args.s-1))
    count_args = _prepare_args_for_counting(
        args.tabix,
        regions_list,
        args.e,
        args.p)
    if args.p > 1:
        with multiprocessing.Pool(args.p) as pool:
            results = pool.starmap(compute_tss_score, count_args)
        results.sort(key=operator.itemgetter("idx"))
    else:
        results = list(itertools.starmap(compute_tss_score, count_args))

    per_barcode_output = f"{args.prefix}.tss_enrichment_barcode_stats.tsv"
    tss_enrichment_plot_fnp = f"{args.prefix}.tss_enrichment_bulk.png"

    with open(per_barcode_output, "w") as out_file:
        # reads_tss_total exists to match the counts produced by ArchR.
        # A read can be counted twice because in two different TSSs.
        # I report the number of reads in TSS without double counting
        # but also the total count like ArchR
        print("barcode\tfragments_promoter\treads_tss\treads_promoter\ttss_enrichment\treads_tss_total", file=out_file)
        for barcode in set([key for result in results for key in result["barcode_statistics"].keys()]):
            merged_signal, merged_statistics = _merge_barcode_results(barcode,
                                                                      [result["barcode_signal"] for result in results],
                                                                      [result["barcode_statistics"] for result in results])
            tss_enrichment, reads_sum = compute_tss_enrichment_barcode(merged_signal)
            print(f"{barcode}\t{merged_statistics[2]}\t{merged_statistics[0]}\t{merged_statistics[1]}\t{tss_enrichment}\t{reads_sum}", file=out_file)

    with open(f"{args.prefix}.tss_score_bulk.txt", "w") as out_file:
        tss_score_bulk = compute_tss_enrichment(_merge_bulk_signal([result["bulk_signal"] for result in results]),
                                                args.window,
                                                tss_enrichment_plot_fnp)
        print(f"tss_enrichment\n{tss_score_bulk}", file=out_file)
