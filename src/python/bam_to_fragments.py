# From Kundaje lab
# https://github.com/kundajelab/ENCODE_scatac/blob/master/workflow/scripts/bam_to_fragments.py

import argparse
import pysam
import sys

def bam_to_frag(in_path, out_path, barcode_tag="CB", shift_plus=4, shift_minus=-4):
    """
    Convert coordinate-sorted BAM file to a fragment file format, while adding Tn5 coordinate adjustment
    BAM should be pre-filtered for PCR duplicates, secondary alignments, and unpaired reads
    Output fragment file is sorted by chr, start, end, barcode
    """

    input = pysam.AlignmentFile(in_path, "rb")
    with open(out_path, "w") as out_file:
        buf = []
        curr_pos = None
        for read in input:
            if read.flag & 16 == 16:
                continue # ignore reverse (coordinate-wise second) read in pair

            chromosome = read.reference_name
            start = read.reference_start + shift_plus
            end = read.reference_start + read.template_length + shift_minus
            cell_barcode = read.get_tag(barcode_tag)
            # assert(read.next_reference_start >= read.reference_start) ####
            data = (chromosome, start, end, cell_barcode, 1)
            pos = (chromosome, start)

            if pos == curr_pos:
                buf.append(data)
            else:
                buf.sort()
                for i in buf:
                    print(*i, sep="\t", file=out_file)
                buf.clear()
                buf.append(data)
                curr_pos = pos

if __name__ == '__main__':

    msg = "Add the description"
    parser = argparse.ArgumentParser(description = msg)

    # Adding optional argument
    parser.add_argument("bam", help = "Path to the coordinate-sorted bam file.")
    parser.add_argument("-o", "--output", help = "Path to the fragments output file.")
    parser.add_argument("--prefix", help = "Prefix for the metrics output file.")
    parser.add_argument("--shift_plus", help = "Tn5 coordinate adjustment for the plus strand.", type = int, default = 4)
    parser.add_argument("--shift_minus", help = "Tn5 coordinate adjustment for the minus strand.", type = int, default = -4)
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
        out_path = f"{prefix}.fragments.tsv"

    bc_tag = args.bc_tag


    bam_to_frag(args.bam, out_path, bc_tag, shift_plus=args.shift_plus, shift_minus=args.shift_minus)
