#!/usr/bin/env python3

"""
This script takes in a bam file, and outputs a txt file containing the number of 
total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode. 
"""

import argparse
import pysam

def parse_arguments():
    parser = argparse.ArgumentParser(description="Get total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode from bam file")
    parser.add_argument("bam_file", help="Filename for input bam file")
    parser.add_argument("bai_file", help="Filename for bam index file")
    parser.add_argument("barcode_metadata_file", help="Filename for output barcode metadata txt file")

    return parser.parse_args()

def add_element(dictionary, key, value):
    if key not in dictionary:
        dictionary[key] = []
    dictionary[key].append(value)

def get_metrics(bam):
    """
    Get barcode metrics from bam file; all counts are only for reads overlapping genes.
    Reported metrics are total counts, UMIs (one UMI counted per unique UMI-gene mapping),
    duplicate counts, percent mitochondrial reads
    """
    total_counts = {}
    genes = {}
    umi_gene = {}
    mitochondrial_counts = {}
    formatted_barcodes = {}
    
    for read in bam.fetch():
        # get barcode; skip read if not present
        barcode = read.get_tag("CB")
        if barcode == "-":
            continue

        # get gene id; skip read if not present
        gene_id = read.get_tag("GX")
        if gene_id == "-":
            continue
        
        # increment total count of barcode    
        total_counts[barcode] = total_counts.get(barcode, 0) + 1
        
        # append gene id to list of genes associated with barcode
        add_element(genes, key=barcode, value=gene_id)
        
        # get umi; add to list of umi-gene combos associated with the barcode
        umi = read.get_tag("UB")
        if umi != "-":
            add_element(umi_gene, key=barcode, value=umi+gene_id)
        
        # if read is mitochondrial, increment mitochondrial count
        if read.reference_name == "chrM":
            mitochondrial_counts[barcode] = mitochondrial_counts.get(barcode, 0) + 1
        
        # format barcode for output (R1,R2,R3,PKR) using barcode from CB tag;
        # read id barcode is not error-corrected so should not be used
        # get PKR from read id rather than workflow "prefix" argmuent to ensure consistency with ATAC
        read_id = read.query_name
        read_id_barcode = read_id.split("_")[1]
        pkr = read_id_barcode.split(",")[3]
        formatted = barcode + "_" + pkr
        formatted_barcodes[barcode] = formatted
    # count unique genes per barcode
    genes_per_barcode = {barcode:len(set(gene_list)) for (barcode, gene_list) in genes.items()}
    
    # count unique umi-gene mappings per barcode
    umis_per_barcode = {barcode:len(set(umi_gene_list)) for (barcode, umi_gene_list) in umi_gene.items()}
    
    # create list with barcodes and associated metrics
    barcode_metadata = []
    for barcode in formatted_barcodes.keys():
        formatted_barcode = formatted_barcodes[barcode]
        total_val = str(total_counts[barcode])
        umi_val = str(umis_per_barcode.get(barcode, 0))
        duplicate_val = str(total_counts[barcode] - umis_per_barcode[barcode])
        gene_val = str(genes_per_barcode[barcode])
        mitochondrial_val = str(round(mitochondrial_counts.get(barcode, 0) / total_counts[barcode] * 100, 2))
        
        metrics = [formatted_barcode, total_val, duplicate_val, umi_val, gene_val, mitochondrial_val]
        
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
    barcode_metadata_file = getattr(args, "barcode_metadata_file")
    
    # load bam file
    bam = pysam.AlignmentFile(bam_file, "rb", index_filename=bai_file)
    
    # get metrics for each barcode
    barcode_metadata = get_metrics(bam)
    
    # write txt file
    write_metadata_file(barcode_metadata, barcode_metadata_file)
    

if __name__ == "__main__":
    main()
