# From Kundaje lab
# https://github.com/kundajelab/ENCODE_scatac/blob/master/workflow/scripts/bam_to_fragments.py

import pysam
import sys

def bam_to_frag(in_path, out_path, shift_plus=4, shift_minus=-4):
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
            cell_barcode = read.get_tag("CB")
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
    try:
        in_path, = sys.argv[1]
        out_path, = sys.argv[2]

        shift_plus = sys.argv[3]
        shift_minus = sys.argv[4]

        bam_to_frag(in_path, out_path, shift_plus=shift_plus, shift_minus=shift_minus)

    except NameError:
        bam_to_frag('/dev/stdin', '/dev/stdout')
