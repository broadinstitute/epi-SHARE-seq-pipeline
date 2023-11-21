# From Kundaje lab
# https://github.com/kundajelab/ENCODE_scatac/blob/master/workflow/scripts/filter_mito.py
#
# Modified by: Eugenio Mattei
# Affiliation: The Broad InstituteOf MIT and Harvard
#
# Changelog:
# 2023/01/20: Now it returns the statistics per barcode
#

import argparse
import pysam
from collections import defaultdict



def filter_mito(in_path, out_path, barcode_tag, cutoff, prefix, threads=1):
    """
    Removes mitochondrial alignments from BAM
    Calculates number of mapped mitochondrial and non-mitochondrial reads (not alignments)
    Assumes mitochondrial chromosome is "chrM"
    """

    infile = pysam.AlignmentFile(in_path, "rb", threads=threads)
    outfile = pysam.AlignmentFile(out_path, "wb", template=infile, threads=threads)
    mito_path = f"{prefix}.mito_only.bam"
    mitofile = pysam.AlignmentFile(mito_path, "wb", template=infile, threads=threads)
    outfile_bulk_metrics = f"{prefix}.mito.bulk-metrics.tsv"
    outfile_barcode_metrics = f"{prefix}.mito.bc-metrics.tsv"

    number_mito = 0
    number_non_mito = 0

    # Initializing the dictionary setting the counts for non-mito and mito.
    barcode_metrics = defaultdict(lambda: [0,0])

    for read in infile.fetch(until_eof=True,multiple_iterators=True):
        if read.reference_name == "chrM":
            if read.flag & 260 == 0: # Alignment is mapped and is primary
                number_mito += 1
                barcode_metrics[read.get_tag(barcode_tag)][1] += 1

        else:
            if read.flag & 260 == 0:
                number_non_mito += 1
                barcode_metrics[read.get_tag(barcode_tag)][0] += 1
            #outfile.write(read)

    # Write the summary metrics
    with open(outfile_bulk_metrics, "w") as fh:
        print("raw_reads_nonmito\traw_reads_mito", file = fh)
        print(f"{number_non_mito}\t{number_mito}", file = fh)

    # Write the metrics per barcode
    with open(outfile_barcode_metrics, "w") as fh:
        # Print header
        print("barcode\traw_reads_nonmito\traw_reads_mito", file = fh)
        for barcode,counts in barcode_metrics.items():
            print(f"{barcode}\t{counts[0]}\t{counts[1]}", file = fh)

    # Write a filtered bam
    for read in infile:
        if read.flag & 260 == 0 and read.reference_name != "chrM" and barcode_metrics[read.get_tag(barcode_tag)][0] > cutoff*2:
            outfile.write(read)
        if read.flag & 260 == 0 and read.reference_name == "chrM":
            mitofile.write(read)

    outfile.close()
    mitofile.close()
    return



if __name__ == '__main__':

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("bam", help = "Path to the coordinate-sorted bam file.")
    parser.add_argument("-o", "--output", help = "Path to the mitochondrial-free bam file.")
    parser.add_argument("-p", help = "Number of threads to use.", type=int, default=1)
    parser.add_argument("--prefix", help = "Prefix for the metrics output file.")
    parser.add_argument("--cutoff", help = "Remove barcodes with a number of fragments less than the cutoff.", type=int, default=1)
    parser.add_argument("--bc_tag", help = "Specify the tag containing the cell barcode.", default="CB")

    # Read arguments from command line
    args = parser.parse_args()

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = args.bam[:-4]

    if args.output:
        out_path = args.output
    else:
        out_path = f"{prefix}.no_mito.bam"

    bc_tag = args.bc_tag

    filter_mito(args.bam, out_path, bc_tag, args.cutoff, prefix, threads=args.p)
