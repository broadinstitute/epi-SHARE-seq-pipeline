#!/usr/bin/env python3

"""
This script takes in a bam file, and outputs a txt file containing the number of
total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode.
"""

import argparse
import logging
import pysam
from collections import defaultdict

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
    total_counts = defaultdict(int)
    genes = defaultdict(set)
    umi_gene = defaultdict(set)
    mitochondrial_counts = defaultdict(int)
    barcodes = set()
    formatted_barcodes = {}

    for read in bam:
        try:
            # get barcode; skip read if not present
            barcode = read.get_tag(barcode_tag)
            if barcode == "-":
                #logging.warning(f"Skipping {read.qname} because the {barcode_tag} tag is empty") slowing down
                continue

            # get gene id; skip read if not present
            gene_id = read.get_tag("GX")
            if gene_id == "-":
                gene_id = read.get_tag("gx")
                if gene_id == "-":
                    #logging.warning(f"Skipping {read.qname} because the GX tag is empty")
                    continue

            # get UMI; skip read if not present
            umi = read.get_tag("UB")
            if umi == "-":
                #logging.warning(f"Skipping {read.qname} because the UB tag is empty")
                continue

            barcodes.add(barcode)

            total_counts[barcode] += 1

            genes[barcode].add(gene_id)

            umi_gene[barcode].add(umi + gene_id)

            if read.reference_name == "chrM":
                mitochondrial_counts[barcode] += 1
        except KeyError:
            logging.error(f"Skipping {read.qname} because one of the tags {barcode_tag},GX, or UB is missing.")

    # count unique genes per barcode
    genes_per_barcode = {barcode:len(gene_set) for (barcode, gene_set) in genes.items()}

    # count unique umi-gene mappings per barcode
    umis_per_barcode = {barcode:len(umi_gene_set) for (barcode, umi_gene_set) in umi_gene.items()}

    # create list with barcodes and associated metrics
    barcode_metadata = []
    for barcode in barcodes:
        total_val = str(total_counts[barcode])
        umi_val = str(umis_per_barcode.get(barcode, 0))
        duplicate_val = str(total_counts[barcode] - umis_per_barcode.get(barcode, 0))
        gene_val = str(genes_per_barcode.get(barcode, 0))
        mitochondrial_val = str(round(mitochondrial_counts.get(barcode, 0) / total_counts[barcode] * 100, 2))
        out_barcode = barcode + "_" + pkr if pkr else barcode

        metrics = [out_barcode, total_val, duplicate_val, umi_val, gene_val, mitochondrial_val]

        barcode_metadata.append(metrics)

    return barcode_metadata

def write_metadata_file(barcode_metadata, output_file):
    fields = ["barcode", "total_counts", "duplicate_counts", "umis", "genes", "percent_mitochondrial"]

    with open(output_file, "w") as f:
        # write header
        f.write("\t".join(fields) + "\n")
        # write rows
        for metrics_list in barcode_metadata:
            f.write("\t".join(metrics_list[:]) + "\n")

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

    # write txt file
    write_metadata_file(barcode_metadata, barcode_metadata_file)

if __name__ == "__main__":

    main()
