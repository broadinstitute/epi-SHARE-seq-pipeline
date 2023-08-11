#!/usr/bin/env python3

"""
This script takes in a bam file, and outputs a txt file containing the number of
total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode.
"""

import argparse
import logging
import numpy as np
import pysam
from collections import defaultdict
from functools import partial

logging.basicConfig(filename='barcode_metadata.log', encoding='utf-8', level=logging.DEBUG)
logging.debug('Creating the barcode metadata for RNA from bam.')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Get total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode from bam file")
    parser.add_argument("bam_file", help="Filename for input bam file")
    parser.add_argument("bai_file", help="Filename for bam index file")
    parser.add_argument("barcode_metadata_file", help="Filename for output barcode metadata txt file")
    parser.add_argument("pkr", help="PKR id for shareseq", default = None, nargs='?')
    parser.add_argument("--barcode_tag", help="PKR id for shareseq", default="CB")

    return parser.parse_args()

def get_metrics(bam, barcode_tag="CB", pkr=None):
    """
    Get barcode metrics from bam file; all counts are only for reads overlapping genes.
    Reported metrics are total counts, UMIs (one UMI counted per unique UMI-gene mapping),
    duplicate counts, genes, percent mitochondrial reads
    """

    # Reads per barcode is a vector of length 2.
    # [x_mito_reads, mito_reads]
    reads_per_barcode = defaultdict(partial(np.zeros,2))
    genes_per_barcode = defaultdict(set)
    mito_genes_per_barcode = defaultdict(set)


    for read in bam:
        try:
            # get barcode; skip read if not present
            barcode = read.get_tag(barcode_tag)
            if barcode == "-":
                #logging.warning(f"Skipping {read.qname} because the {barcode_tag} tag is empty") slowing down
                continue
                
            # get UMI; skip read if not present
            umi = read.get_tag("UB")
            if umi == "-":
                #logging.warning(f"Skipping {read.qname} because the UB tag is empty")
                continue

            reads_per_barcode[barcode][0] += 1
            
            if read.reference_name == "chrM":
                reads_per_barcode[barcode][1] += 1

        except KeyError:
            #logging.error(f"Skipping {read.qname} because one of the tags {barcode_tag},GX, or UB is missing.")
            continue

        try:
            # get gene id; skip read if not present
            gene_id = read.get_tag("GX")
            if gene_id == "-":
                #logging.warning(f"Skipping {read.qname} because the GX tag is empty")
                continue

            if read.reference_name == "chrM":
                mito_genes_per_barcode[barcode].add(gene_id)
            else:
                genes_per_barcode[barcode].add(gene_id)
        except KeyError:
            #logging.error(f"Skipping {read.qname} because one of the tags {barcode_tag},GX, or UB is missing.")
            continue

    # create list with barcodes and associated metrics
    barcode_metadata = []
    for barcode, reads_vector in reads_per_barcode.items():
        # Reminder that reads_vector is [x_mito_reads, mito_reads].
        genes = len(genes_per_barcode[barcode])
        mito_genes = len(mito_genes_per_barcode[barcode])
        fraction_mitochondrial_reads = round(reads_vector[1]/np.sum(reads_vector) * 100, 2)
        out_barcode = barcode + "_" + pkr if pkr else barcode

        metrics = list(map(str,[out_barcode, int(np.sum(reads_vector)), reads_vector[0], reads_vector[1], genes+mito_genes, genes, mito_genes, fraction_mitochondrial_reads]))

        barcode_metadata.append(metrics)

    return (barcode_metadata)


def write_metadata_file(barcode_metadata, output_file):
    fields = ["barcode", "total_counts", "reads_non_mito", "reads_mito", "genes", "genes_non_mito", "genes_mito", "percent_mitochondrial"]

    with open(output_file, "w") as f:
        # write header
        f.write("\t".join(fields) + "\n")
        # write rows
        for metrics_list in barcode_metadata:
            f.write("\t".join(metrics_list) + "\n")

def main():
    # get arguments
    args = parse_arguments()
    bam_file = getattr(args, "bam_file")
    bai_file = getattr(args, "bai_file")
    pkr = getattr(args, "pkr")
    barcode_tag = getattr(args, "barcode_tag")
    barcode_metadata_file = getattr(args, "barcode_metadata_file")

    # load bam file
    bam = pysam.AlignmentFile(bam_file, "rb", index_filename=bai_file)

    # get metrics for each barcode
    barcode_metadata = get_metrics(bam, barcode_tag, pkr)

    # write metadata file
    write_metadata_file(barcode_metadata, barcode_metadata_file)
    
if __name__ == "__main__":

    main()
