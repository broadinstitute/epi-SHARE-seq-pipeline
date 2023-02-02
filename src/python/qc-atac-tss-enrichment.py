#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University
# modified from plotV_vC.py (Alicia)

# Will make a V-plot from bed regions

import argparse
import os
import pysam
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
                                center=True,
                                shift_plus=4,
                                shift_minus=-4,
                                barcode_tag = "CB",
                                mapq_threshold = 30):
    """
    The function that does the counting for each input BAM/SAM file.
    Fixme: there are some redundant parameters here.. feature_type, id_attribute, additional_attributes
    Parameters
    ----------
    bam_filename : str
        Path to the SAM/BAM file containing the mapped reads.
        File needs to be coordinate-sorted and indexed.
    tss_list : array
        TODO check the type of this parameter.
        Supplied by HTSeq.make_feature_genomicarrayofsets
    upstream : array
        TODO check the type of this parameter.
        Supplied by HTSeq.make_feature_genomicarrayofsets
    downstream : str
        Can only be either 'pos' or 'name'. Sorting order of <alignment_file>.
    mapq_threshold : int
        The number of reads allowed to stay in memory until mates are found.
        Used when <alignment_file> is paired end sorted by position.
    center : int
        The number of reads allowed to stay in memory until mates are found.
        Used when <alignment_file> is paired end sorted by position.
    shift_plus : int
        Whether the data to be aligned is from a strand-specific assay.
        Option is yes, no, reverse.
        reverse means yes with reversed strand interpretation.
    shift_minus : str
        Mode to handle reads overlapping more than one feature.
        Choices: union, intersection-strict, intersection-nonempty.
    barcode_tag : str
        Whether and how to score reads that are not uniquely aligned or
        ambiguously assigned to features.
        Choices: none, all, fraction, random.
    Returns
    -------
    Dictionary
        Key: Barcode
        Value: array[TSS-upstream,TSS+downstream] containing the number of fragments per basepair.
    """
    promoter_size = upstream+downstream
    counts_dict = defaultdict(lambda: np.zeros(300))
    counts_dict_bulk = np.zeros(promoter_size)
    fragments_in_promoter_counter = defaultdict(int)

    bamfile = pysam.Samfile(bam_filename, "rb")

    for tss in tss_list:
        # tss example: ["chr", "start", "end", "strand"]
        # Find the center of the TSS
        tss_center = int(tss[1])+(int(tss[2])-int(tss[1]))/2
        # Create the promoter region by adding the upstream and downstream.
        promoter_start = tss_center - upstream
        promoter_end = tss_center + downstream
        promoter_strand = tss[3]
        # Find all the fragments overlapping the promoter.
        # Why is this adding an extra region?
        for read in bamfile.fetch(str(tss[0]), max(0,promoter_start), promoter_end):
            #check mapping quality
            if read.mapq < mapq_threshold or read.flag & 16 == 16 or read.get_tag(barcode_tag) !="TCATCCTAGCAGAGCCTCGAGCGT":
                continue # Ignore low quality reads and reverse (coordinate-wise second) read in pair.

            fragment_start = read.pos +  shift_plus
            fragment_length = abs(read.tlen) - shift_plus #? - shift_plus + shift_minus
            fragment_end = fragment_start + fragment_length
            fragment_center = fragment_start + (fragment_length/2)

            # Increment the counter for the specific barcode.
            barcode = read.get_tag(barcode_tag)
            fragments_in_promoter_counter[barcode] += 1

            # Update the array with the counts around the promoter.
            if center:
                # The fragment spans two bp so we need to update them both and give each half the weight
                if fragment_length %2 == 0:
                    _add_read_to_dictionary(counts_dict_bulk, counts_dict, fragment_center, promoter_start, promoter_end, promoter_strand, barcode, 1)
                else:
                    _add_read_to_dictionary(counts_dict_bulk, counts_dict, fragment_center, promoter_start, promoter_end, promoter_strand, barcode, 0.5)
                    _add_read_to_dictionary(counts_dict_bulk, counts_dict, fragment_center+1, promoter_start, promoter_end, promoter_strand, barcode, 0.5)
            else:
                _add_read_to_dictionary(counts_dict_bulk, counts_dict, fragment_start, promoter_start, promoter_end, promoter_strand, barcode, 1)
                _add_read_to_dictionary(counts_dict_bulk, counts_dict, fragment_end, promoter_start, promoter_end, promoter_strand, barcode, 1)

    return counts_dict_bulk, counts_dict, fragments_in_promoter_counter

def _add_read_to_dictionary(fragment_counter,
                            fragment_counter_barcodes,
                            fragment_position,
                            promoter_start,
                            promoter_end,
                            promoter_strand,
                            barcode,
                            weight):
    if  fragment_position >= promoter_start and fragment_position <= promoter_end-1:
        # We are adding to an array that covers each basepair from promoter_start to promoter_end.
        # This converts the genomic coordinates to the index that we need to update in the array.
        index = int(fragment_position - promoter_start)
        max_window = len(fragment_counter)

        add_count_to_barcode = False

        if index >= 0 and index <= 99:
            index_reduced = index
            add_count_to_barcode = True
        if index >= max_window-100:
            index_reduced = index - (max_window-300)
            add_count_to_barcode = True
        if index >= (max_window/2)-50 and index <= (max_window/2)+50:
            # Find the center of the region with max_window/2+50
            index_reduced = int(index + 100 - ((max_window/2)-50))
            add_count_to_barcode = True

        # If the TSS is in the negative strand we need to add from the end.
        if promoter_strand == "-":
            index_reverse = max_window - index - 1
            fragment_counter[index_reverse] += weight

            if add_count_to_barcode:
                index_reduced_reversed = 300 - index_reduced - 1
                fragment_counter_barcodes[barcode][index_reduced_reversed] += weight
        else:
            fragment_counter[index] += weight

            if add_count_to_barcode:
                fragment_counter_barcodes[barcode][index_reduced] += weight

    return

def compute_tss_enrichment(array_counts, window_size):

    size = len(array_counts)
    # We want to normalize the plot using the signal in the first and last 200 bps.
    #normalization_factor = np.mean(array_counts[np.r_[1:100, size-100:size]])
    normalization_factor = np.mean(array_counts[np.r_[0:100, size-100:size]])
    raw_signal = array_counts/normalization_factor
    # Smooth using a window
    smoothed_signal = np.convolve(array_counts, np.ones(window_size),'same')/window_size/normalization_factor
    fig=plt.figure(figsize=(8.0, 5.0))
    plt.plot(raw_signal,'k.')
    plt.plot(smoothed_signal,'r')
    plt.xlabel('Position relative to center')
    plt.ylabel('Insertions')
    fig.savefig('test-flag-correct-promoter-1000.png')
    plt.close(fig)
    return max(raw_signal)

def compute_tss_enrichment_cell(array_counts):
    size = len(array_counts)
    # We want to normalize the plot using the signal in the first and last 200 bps.
    #normalization_factor = np.mean(array_counts[np.r_[1:100, size-100:size]])
    normalization_factor = np.mean(array_counts[np.r_[0:100, 200:300]])
    return np.mean(array_counts[100:200])/normalization_factor


if __name__ == '__main__':

    #args = _parse_sanitize_cmdline_arguments()

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("bam", help= "Path to the coordinate-sorted bam file.")
    #parser.add_argument("-o", "--output", help= "Path to the mitochondrial-free bam file.")
    parser.add_argument("-e", help= "Number of bases to extend to each side. (default= 1000)", type= int, default= 2000)
    parser.add_argument("--ends", help= "Use the ends of the fragment to count. Default is to use the center.", action= "store_true")
    parser.add_argument("--mapq_threshold", help= "Filter reads with a mapq value lower than the threshold.", type= int, default= 30)
    parser.add_argument("-s", help="Column with strand information; 1-based. (default= 4)", type= int, default= 4)
    parser.add_argument("-c", help= "Number of cpus to use. (default= max/2)")
    parser.add_argument("-w", "--window", help= "Smoothing window size for plotting. (default= 20)", type= int, default= 20)
    parser.add_argument("--bc_tag", help = "Specify the tag containing the cell barcode.", default="CB")
    parser.add_argument("--prefix", help = "Prefix for the metrics output file.")
    parser.add_argument("--tss", help= "TSS bed file")
    parser.add_argument("--fragment_cutoff", help= "TSS bed file", type= int, default= 10)


    # Read arguments from command line
    args = parser.parse_args()

    #cpus = len(os.sched_getaffinity(0))/2

    tss_list = np.loadtxt(args.tss, 'str', usecols = (0,1,2,args.s-1))

    bulk_counts, barcode_counts, stats = count_fragments_in_promoter(args.bam,
                                tss_list,
                                flank= args.e,
                                center=False,
                                shift_plus=4,
                                shift_minus=-4,
                                barcode_tag = args.bc_tag,
                                mapq_threshold = args.mapq_threshold
                                )

    for barcode, fragments_in_promoter in stats.items():
        if fragments_in_promoter >= args.fragment_cutoff:
            score = compute_tss_enrichment_cell(barcode_counts[barcode])
            print(barcode_counts[barcode])
            print(f"{barcode} {score}")
    print(bulk_counts[0:100])
    print(bulk_counts[3900:4000])
    print(bulk_counts[1950:2050])
    print(np.mean(bulk_counts[1950:2050]))
    print(np.mean(bulk_counts[1950:2050])/0.2)
    print(f"bulk_score:{compute_tss_enrichment(bulk_counts, args.window)}")

    # It takes 30 seconds to run on a small sample. I don't think it needs to be parallel.
    # split bedfile into chunks
    #tss_number = len(tss_list)
    #chunksize = tss_number/int(cpus)
    #starts = range(0,tss_number, chunksize)
    #if cpus > 1:
    #    with multiprocessing.Pool(cpus) as pool:
    #        results = pool.starmap(count_reads_single_file, count_args)
    #else:
    #    results = list(itertools.starmap(count_fragments_in_promoter, count_args))


    #if args.prefix:
    #    prefix = args.prefix
    #else:
    #    prefix = args.bam[:-4]

    #if args.output:
    #    out_path = args.output
    #else:
    #    out_path = f"{prefix}.no_mito.bam"

    #bc_tag = args.bc_tag
