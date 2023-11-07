#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make uncorrected R1 and R2 FASTQs from unmapped BAM file.
"""

import argparse
import pysam
from utils import check_putative_barcode, create_barcode_dicts


def parse_arguments():
    parser = argparse.ArgumentParser(description="Make R1 and R2 FASTQs from unmapped BAM file")
    parser.add_argument("bam_file", help="Filename for unmapped BAM file")
    parser.add_argument("pkr", help="PKR name")
    parser.add_argument("prefix", help="Prefix for naming output FASTQ files")
    parser.add_argument("r1_barcode_set_file", help="File containing biosample splits in R1 barcodes, one split per lane")
    parser.add_argument("--r2_barcode_file", help="File containing R2 barcodes, one line")
    parser.add_argument("--r3_barcode_file", help="File containing R3 barcodes, one line")

    return parser.parse_args()


# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}


def reverse_complement(sequence):
    """
    Return reverse complement of DNA sequence.
    """
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])


def write_read(fastq, read):
    """
    Write read to open FASTQ file.
    """
    info = {'index': int(not read.is_read1) + 1,
            'name':  read.query_name}

    if read.is_reverse:
        info.update({'quality':  read.qual[::-1],
                     'sequence': reverse_complement(read.query_sequence)})
    else:
        info.update({'quality':  read.qual,
                     'sequence': read.query_sequence})

    fastq.write('@{name}\n{sequence}\n+\n{quality}\n'.format(**info))


def create_barcode_subset_dict(file_path):
    """
    Create dictionary mapping barcodes to their R1 barcode subsets
    """
    with open(file_path) as f:
        barcode_subset_dict = dict()  # {barcode: subset}
        for line in f:
            line = line.strip().split()
            subset = line[0]  # first entry of line is subset name
            for barcode in line[1:]:
                barcode_subset_dict[barcode] = subset
    return barcode_subset_dict


def write_fastqs(bam_file,
                 read_1_pointers,
                 read_2_pointers,
                 r1_barcode_subset_dict,
                 r1_barcode_exact_dict,
                 r1_barcode_mismatch_dict,
                 prefix):
    """
    Get reads from BAM file and error-correct R1 barcode to determine which R1 barcode subset's
    FASTQ to write the read to. Barcodes written into FASTQs are not error-corrected.

    read_1_pointers is a dictionary mapping R1 barcode subset names to corresponding
    R1 FASTQ filenames to write to.
    read_2_pointers is a dictionary mapping R1 barcode subset names to corresponding
    R2 FASTQ filenames to write to.
    r1_barcode_matches is a dictionary mapping R1 barcodes to their potential matches
    (exact and 1bp mismatch)

    Outputs raw R1 and R2 FASTQs, as well as QC txt file containing R1 barcode statistics.
    """
    bam = pysam.Samfile(bam_file, "rb", check_sq=False)
    query_name = read_1 = read_2 = None
    exact_match = nonexact_match = nonmatch = poly_g_barcode = 0

    for read in bam:
        if read.is_read1:
            # save and continue processing
            read_1 = read
            query_name = read_1.query_name

        else:
            # check that query names are the same
            if read.query_name == query_name:
                read_2 = read
                barcode_tag = read.get_tag("RX")
                quality_tag = read.get_tag("QX")

                # append barcode and corresponding quality to read 2
                read_2_quality = read_2.qual  # save to variable; defaults to NoneType once sequence is modified
                read_2.query_sequence += barcode_tag
                read_2.qual = read_2_quality + quality_tag

                # get 10bp sequence containing R1 barcode (additional 1bp padding used for checking shifts)
                r1_barcode_window = barcode_tag[14:24]
                # get error-corrected R1 barcode
                r1_barcode, quality, exact = check_putative_barcode(r1_barcode_window,
                                                                    quality_tag,
                                                                    r1_barcode_exact_dict,
                                                                    r1_barcode_mismatch_dict)

                # write reads to appropriate FASTQ files by checking which R1 barcode subset the corrected R1 barcode belongs to
                if r1_barcode in r1_barcode_subset_dict.keys():
                    write_read(read_1_pointers[r1_barcode_subset_dict[r1_barcode]], read_1)
                    write_read(read_2_pointers[r1_barcode_subset_dict[r1_barcode]], read_2)

                # increment QC counters
                if r1_barcode:
                    if exact:
                        exact_match += 1
                    else:
                        nonexact_match += 1
                elif "G"*8 in r1_barcode_window:
                    poly_g_barcode += 1
                else:
                    nonmatch += 1

    # write QC stats
    with open(f"{prefix}_R1_barcode_qc.txt", "w") as f:
        f.write("%s\t%s\t%s\t%s\t%s\n" % (prefix, exact_match, nonexact_match, nonmatch, poly_g_barcode))


def main():
    args = parse_arguments()
    bam_file = getattr(args, "bam_file")
    pkr = getattr(args, "pkr")
    prefix = getattr(args, "prefix")
    r1_barcode_set_file = getattr(args, "r1_barcode_set_file")
    r2_barcode_file = getattr(args, "r2_barcode_file")
    r3_barcode_file = getattr(args, "r3_barcode_file")

    # create dictionary mapping R1 barcodes to their R1 barcode subsets
    r1_barcode_subset_dict = create_barcode_subset_dict(r1_barcode_set_file)
    # create dictionaries of R1 barcode exact matches and 1bp mismatches
    r1_barcode_exact_dict, r1_barcode_mismatch_dict = create_barcode_dicts(r1_barcode_subset_dict.keys())

    # read R2 and R3 barcode files; if not passed in, use R1 barcodes
    if r2_barcode_file:
        r2_barcode_subset_dict = create_barcode_subset_dict(r2_barcode_file)
        r2_barcodes = r2_barcode_subset_dict.keys()
    else:
        r2_barcodes = r1_barcode_subset_dict.keys()

    if r3_barcode_file:
        r3_barcode_subset_dict = create_barcode_subset_dict(r2_barcode_file)
        r3_barcodes = r3_barcode_subset_dict.keys()
    else:
        r3_barcodes = r1_barcode_subset_dict.keys()

    # reads are written into one R1 FASTQ and one R2 FASTQ per R1 barcode subset;
    # create dictionaries of file pointers associated with each R1 barcode subset,
    # make whitelist for each R1 barcode subset
    read_1_pointers = dict()
    read_2_pointers = dict()
    for subset in set(r1_barcode_subset_dict.values()):
        fp = open(prefix + "_" + pkr + "_" + subset + "_R1.fastq", "w")
        read_1_pointers[subset] = fp
        fp = open(prefix + "_" + pkr + "_" + subset + "_R2.fastq", "w")
        read_2_pointers[subset] = fp
        # get possible combinations of R1R2R3, write to whitelist
        r1_barcodes = [k for k, v in r1_barcode_subset_dict.items() if v == subset]
        whitelist_barcodes = [r1+r2+r3 for r1 in r1_barcodes for r2 in r2_barcodes for r3 in r3_barcodes]
        with open(prefix + "_" + pkr + "_" + subset + "_whitelist.txt", "w") as f:
            f.write("\n".join(whitelist_barcodes))

    # write reads to FASTQs
    write_fastqs(bam_file,
                 read_1_pointers,
                 read_2_pointers,
                 r1_barcode_subset_dict,
                 r1_barcode_exact_dict,
                 r1_barcode_mismatch_dict,
                 prefix)


if __name__ == "__main__":
    main()
