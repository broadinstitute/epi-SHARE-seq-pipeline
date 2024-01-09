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
    parser.add_argument("subpool", help="Cellular subpool name")
    parser.add_argument("prefix", help="Prefix for naming output FASTQ files")
    parser.add_argument("r1_barcode_set_file", help="File containing round 1 barcodes of plate subsets, one subset per line")
    parser.add_argument("--r2_barcode_file", help="File containing round 2 barcodes, one line")
    parser.add_argument("--r3_barcode_file", help="File containing round 3 barcodes, one line")

    return parser.parse_args()


# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}


def reverse_complement(sequence):
    """Return reverse complement of DNA sequence."""
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])


def write_read(fastq, read):
    """Write read to open FASTQ file."""
    info = {"index": int(not read.is_read1) + 1,
            "name": read.query_name}

    if read.is_reverse:
        info.update({"quality": read.qual[::-1],
                     "sequence": reverse_complement(read.query_sequence)})
    else:
        info.update({"quality":  read.qual,
                     "sequence": read.query_sequence})

    fastq.write("@{name}\n{sequence}\n+\n{quality}\n".format(**info))


def create_barcode_subset_dict(file_path):
    """Create dictionary mapping barcodes to their R1 barcode subsets."""
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
                 r1_barcode_mismatch_dict):
    """
    Write BAM reads into R1 and R2 FASTQs based on corrected round 1 barcode.

    Get reads from unaligned BAM file and error-correct round 1 barcode to
    determine which round 1 barcode subset FASTQs to write the read to.
    Writes an R1 and R2 FASTQ for each round 1 subset. Barcodes written into
    FASTQs are ~NOT~ error-corrected.

    Inputs:
    read_1_pointers: dictionary mapping round 1 barcode subset names to
    corresponding R1 FASTQ filenames.
    read_2_pointers: dictionary mapping round 1 barcode subset names to
    corresponding R2 FASTQ filenames.
    r1_barcode_subset_dict: dictionary mapping round 1 barcodes to their
    round 1 barcode subsets.
    r1_barcode_exact_dict: dictionary mapping round 1 barcodes to their exact
    matches.
    r1_barcode_mismatch_dict: dictionary mapping round 1 barcodes to their
    one-base mismatches.

    Outputs:
    library_qcs: list containing total number of reads, number of reads with
    accepted round 1 barcodes, number of reads discarded due to poly-G round 1
    barcodes, and number of reads discarded due to unaccepted-but-not-poly-G
    round 1 barcodes.
    r1_barcode_qcs: dictionary mapping round 1 barcodes to the total number of
    times that barcode was observed, as well as how many times that barcode
    was recovered through exact matching, left shifting, right shifting, and
    one-base mismatching.
    """
    # initialize dict of per-barcode QC metrics;
    # [exact match, left shift, right shift, mismatch]
    r1_barcode_qcs = {barcode: [0, 0, 0, 0] for barcode in r1_barcode_subset_dict.keys()}
    # initialize library QC metrics
    r1_match = r1_poly_g = r1_nonmatch = 0

    bam = pysam.Samfile(bam_file, "rb", check_sq=False)
    query_name = read_1 = read_2 = None

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

                # get 10bp sequence containing round 1 barcode (additional 1bp padding used for checking shifts)
                r1_barcode_window = barcode_tag[14:24]
                # get error-corrected round 1 barcode
                r1_barcode, quality, correction_type = check_putative_barcode(r1_barcode_window,
                                                                              quality_tag,
                                                                              r1_barcode_exact_dict,
                                                                              r1_barcode_mismatch_dict)

                if r1_barcode:
                    # write reads to appropriate FASTQ files by checking which
                    # round 1 barcode subset the corrected round 1 barcode belongs to
                    write_read(read_1_pointers[r1_barcode_subset_dict[r1_barcode]], read_1)
                    write_read(read_2_pointers[r1_barcode_subset_dict[r1_barcode]], read_2)

                    r1_match += 1

                    if correction_type == "E":
                        r1_barcode_qcs[r1_barcode][0] += 1
                    elif correction_type == "L":
                        r1_barcode_qcs[r1_barcode][1] += 1
                    elif correction_type == "R":
                        r1_barcode_qcs[r1_barcode][2] += 1
                    elif correction_type == "M":
                        r1_barcode_qcs[r1_barcode][3] += 1

                elif "GGGGGGGG" in r1_barcode_window:
                    r1_poly_g += 1

                else:
                    r1_nonmatch += 1

    total_reads = r1_match + r1_poly_g + r1_nonmatch
    library_qcs = [total_reads, r1_match, r1_poly_g, r1_nonmatch]
    r1_barcode_qcs = {k: v.append(sum(v)) for k, v in r1_barcode_qcs.items()}  # append total

    return library_qcs, r1_barcode_qcs


def write_qc_stats(prefix, subpool, library_qcs, r1_barcode_qcs):
    """Write out library and per-barcode QC statistics"""
    library_qc_fields = ["subpool", "library", "total_reads", "r1_match", "r1_poly_G", "r1_nonmatch"]
    with open(f"{prefix}_{subpool}_library_qc.txt", "w") as f:
        f.write("\t".join(library_qc_fields) + "\n")
        f.write(f"{subpool}\t{prefix}\t" + "\t".join(map(str, library_qcs)))

    r1_barcode_qc_fields = ["r1_barcode", "exact_match", "left_shift", "right_shift", "mismatch", "total"]
    with open(f"{prefix}_{subpool}_R1_barcode_qc.txt", "w") as f:
        f.write("\t".join(r1_barcode_qc_fields) + "\n")
        for k, v in r1_barcode_qcs.items():
            f.write(f"{k}\t" + "\t".join(map(str, v)))


def main():
    args = parse_arguments()
    bam_file = getattr(args, "bam_file")
    subpool = getattr(args, "subpool")
    prefix = getattr(args, "prefix")
    r1_barcode_set_file = getattr(args, "r1_barcode_set_file")
    r2_barcode_file = getattr(args, "r2_barcode_file")
    r3_barcode_file = getattr(args, "r3_barcode_file")

    # create dictionary mapping round 1 barcodes to their round 1 barcode subsets
    r1_barcode_subset_dict = create_barcode_subset_dict(r1_barcode_set_file)
    # create dictionaries of round 1 barcode exact matches and 1bp mismatches
    r1_barcode_exact_dict, r1_barcode_mismatch_dict = create_barcode_dicts(r1_barcode_subset_dict.keys())

    # read round 2 and round 3 barcode files; if not passed in, use round 1 barcodes
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

    # reads are written into one R1 FASTQ and one R2 FASTQ per round 1 barcode subset;
    # create dictionaries of file pointers associated with each round 1 barcode subset,
    # make whitelist for each round 1 barcode subset
    read_1_pointers = dict()
    read_2_pointers = dict()
    for subset in set(r1_barcode_subset_dict.values()):
        fp = open(f"{prefix}_{subpool}_{subset}_R1.fastq", "w")
        read_1_pointers[subset] = fp
        fp = open(f"{prefix}_{subpool}_{subset}_R2.fastq", "w")
        read_2_pointers[subset] = fp

        # get possible combinations of R1R2R3, write to whitelist
        r1_barcodes = [k for k, v in r1_barcode_subset_dict.items() if v == subset]
        whitelist_barcodes = [r1+r2+r3 for r1 in r1_barcodes for r2 in r2_barcodes for r3 in r3_barcodes]
        with open(f"{prefix}_{subpool}_{subset}_whitelist.txt", "w") as f:
            f.write("\n".join(whitelist_barcodes))

    # write reads to FASTQs, get QC statistics
    library_qcs, r1_barcode_qcs = write_fastqs(bam_file,
                                               read_1_pointers,
                                               read_2_pointers,
                                               r1_barcode_subset_dict,
                                               r1_barcode_exact_dict,
                                               r1_barcode_mismatch_dict)

    write_qc_stats(prefix, subpool, library_qcs, r1_barcode_qcs)


if __name__ == "__main__":
    main()
