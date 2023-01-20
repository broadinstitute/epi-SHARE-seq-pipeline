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



def filter_mito(in_path, out_path, barcode_tag, prefix):
    """
    Removes mitochondrial alignments from BAM
    Calculates number of mapped mitochondrial and non-mitochondrial reads (not alignments)
    Assumes mitochondrial chromosome is "chrM"
    """

    infile = pysam.AlignmentFile(in_path, "rb")
    outfile = pysam.AlignmentFile(out_path, "wb", template=infile)
    outfile_bulk_metrics = f"{prefix}.bulk-metrics.tsv"
    outfile_barcode_metrics = f"{prefix}.bc-metrics.tsv"

    number_mito = 0
    number_non_mito = 0

    # Initializing the dictionary setting the counts for non-mito and mito.
    barcode_metrics = defaultdict(lambda: [0,0])

    for read in infile:
        if read.reference_name == "chrM":
            if read.flag & 260 == 0: # Alignment is mapped and is primary
                number_mito += 1
                barcode_metrics[read.get_tag("CB")][1] += 1

        else:
            if read.flag & 260 == 0:
                number_non_mito += 1
                barcode_metrics[read.get_tag("CB")][0] += 1
            outfile.write(read)


    # Write the summary metrics
    with open(outfile_bulk_metrics, "w") as fh:
        print("Non-Mitochondrial\tMitochondrial", file = fh)
        print(f"{number_non_mito}\t{number_mito}", file = fh)

    # Write the metrics per barcode
    with open(outfile_barcode_metrics, "w") as fh:
        # Print header
        print("Barcode\tNon-Mitochondrial\tMitochondrial", file = fh)
        for barcode,counts in barcode_metrics.items():
            print(f"{barcode}\t{counts[0]}\t{counts[1]}", file = fh)



if __name__ == '__main__':

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("bam", help = "Path to the coordinate-sorted bam file.")
    parser.add_argument("-o", "--output", help = "Path to the mitochondrial-free bam file.")
    parser.add_argument("--prefix", help = "Prefix for the metrics output file.")
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

    filter_mito(args.bam, out_path, bc_tag, prefix)
