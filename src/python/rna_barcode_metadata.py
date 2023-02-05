#!/usr/bin/env python3

"""
This script takes in a bam file, and outputs a txt file containing the number of 
total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode. 
"""

import argparse
from collections import defaultdict
import pysam

def parse_arguments():
    parser = argparse.ArgumentParser(description="Get total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode from bam file")
    parser.add_argument("chemistry", choices=["shareseq", "10x_v2", "10x_v3", "splitseq"], help="Method chemistry")
    parser.add_argument("bam_file", help="Filename for input bam file")
    parser.add_argument("bai_file", help="Filename for bam index file")
    parser.add_argument("barcode_metadata_file", help="Filename for output barcode metadata txt file")

    return parser.parse_args()

def get_metrics(bam, chemistry):
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
    
    for read in bam.fetch():
        # get barcode; skip read if not present
        barcode = read.get_tag("CB")
        if barcode == "-":
            continue
        
        # get gene id; skip read if not present
        gene_id = read.get_tag("GX")
        if gene_id == "-":
            continue
        
        # get UMI; skip read if not present
        umi = read.get_tag("UB")
        if umi == "-":
            continue
        
        total_counts[barcode] += 1
        
        genes[barcode].add(gene_id)
        
        umi_gene[barcode].add(umi + gene_id)
        
        if read.reference_name == "chrM":
            mitochondrial_counts[barcode] += 1
        
        # SHARE-seq: format barcode for output (R1,R2,R3,PKR) using barcode from CB tag;
        # read id barcode is not error-corrected so should not be used
        # get PKR from read id rather than workflow "prefix" argmuent to ensure consistency with ATAC
        if chemistry == "shareseq":
            read_id = read.query_name
            read_id_barcode = read_id.split("_")[1]
            pkr = read_id_barcode.split(",")[3]
            formatted_barcode = barcode[:8] + "," + barcode[8:16] + "," + barcode[16:] + "," + pkr
            barcodes.add(formatted_barcode)
        
        # all other chemistries: no additional barcode formatting
        else:
            barcodes.add(barcode)
        
    # count unique genes per barcode
    genes_per_barcode = {barcode:len(gene_set) for (barcode, gene_set) in genes.items()}
    
    # count unique umi-gene mappings per barcode
    umis_per_barcode = {barcode:len(umi_gene_set) for (barcode, umi_gene_set) in umi_gene.items()}
    
    # create list with barcodes and associated metrics
    barcode_metadata = []
    for barcode in barcodes:
        total_val = str(total_counts[barcode])
        umi_val = str(umis_per_barcode.get(barcode, 0))
        duplicate_val = str(total_counts[barcode] - umis_per_barcode[barcode])
        gene_val = str(genes_per_barcode[barcode])
        mitochondrial_val = str(round(mitochondrial_counts.get(barcode, 0) / total_counts[barcode] * 100, 2))
        
        metrics = [barcode, total_val, duplicate_val, umi_val, gene_val, mitochondrial_val]
        
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
    chemistry = getattr(args, "chemistry")
    bam_file = getattr(args, "bam_file")
    bai_file = getattr(args, "bai_file")
    barcode_metadata_file = getattr(args, "barcode_metadata_file")
    
    # load bam file
    bam = pysam.AlignmentFile(bam_file, "rb", index_filename=bai_file)
    
    # get metrics for each barcode
    barcode_metadata = get_metrics(bam, chemistry)
    
    # write txt file
    write_metadata_file(barcode_metadata, barcode_metadata_file)
    
if __name__ == "__main__":
    main()
