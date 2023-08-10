

# Count the number of unique and duplicate fragments per barcode using the DS and DI tag from picard

import argparse
import pysam
import sys

from collections import defaultdict

def count_duplicates(in_path, out_path, barcode_tag="CB"):
    """
    """
    # Dictionary holding the unique and duplicate count per barcode
    counter = defaultdict(lambda: [0,0])
    input = pysam.AlignmentFile(in_path, "rb")
    for read in input:
        cell_barcode = read.get_tag(barcode_tag)
        if read.flag & 1024 == 1024:
            counter[cell_barcode][1] += 1
        else:
            counter[cell_barcode][0] += 1

    with open(out_path, "w") as out_file:
        print("barcode\treads_unique\treads_duplicate\tpct_duplicates", file=out_file)
        for barcode, counts_vector in counter.items():
            print(f"{barcode}\t{counts_vector[0]}\t{counts_vector[1]}\t{round(counts_vector[1]/(counts_vector[0]+counts_vector[1])*100,1)}", file=out_file)

if __name__ == '__main__':

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("bam_wdup", help = "Path to the coordinate-sorted bam file with duplicates marked but nor removed.")
    parser.add_argument("-o", "--output", help = "Path to the fragments output file.")
    parser.add_argument("--prefix", help = "Prefix for the metrics output file.")
    parser.add_argument("--bc_tag", help = "Specify the tag containing the cell barcode.", default="CB")

    # Read arguments from command line
    args = parser.parse_args()

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = args.bam_wdup[:-4]

    if args.output:
        out_path = args.output
    else:
        out_path = f"{prefix}.duplicate.stats.tsv"

    bc_tag = args.bc_tag


    count_duplicates(args.bam_wdup, out_path, bc_tag)
