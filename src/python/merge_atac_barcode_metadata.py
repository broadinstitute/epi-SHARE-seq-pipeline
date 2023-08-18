#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge ATAC barcode metadata TSV files
"""

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge ATAC barcode metadata TSV files")
    parser.add_argument("output_file", help="Filename for output merged barcode metadata TSV file")
    parser.add_argument("metadata_files", nargs="*", help="Filenames for barcode metadata TSV files to be merged")
    return parser.parse_args()

def get_split_lines(file_name, delimiter, skip=0):
    """
    Read file contents and yield generator with line entries
    """
    with open(file_name, "rt") as f:
        for i in range(skip):
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)
            
def merge_barcode_metadata(metadata_files):
    """
    Read files and store metadata for each barcode.
    Sum values for duplicate barcodes. 
    Compute pct_duplicates, pct_reads_promoter, pct_reads_peaks, pct_mito_reads.
    """
    metadata_dict = {}
    for metadata_file in metadata_files:
        metadata = get_split_lines(metadata_file, delimiter="\t")
        next(metadata)
        for line in metadata:
            if line[0] in metadata_dict.keys(): 
                # sum if barcode already in dict
                old_list = metadata_dict[line[0]]
                new_list = [int(float(x)) for x in line[1:]] # can't convert '0.5' to int; need to convert to float first
                metadata_dict[line[0]] = [sum(x) for x in zip(old_list, new_list)]
            else:
                metadata_dict[line[0]] = [int(float(x)) for x in line[1:]]
    
    # compute percentages
    for k,v in metadata_dict.items():
        pct_duplicates = round(v[6]/(v[5]+v[6])*100, 1)
        pct_reads_promoter = round(v[2]/v[5]*100, 1)
        pct_reads_peaks = round(v[8]/v[5]*100, 1)
        pct_mito_reads = round(v[11]/(v[10]+v[11])*100, 1)
        v[7] = pct_duplicates
        v[12] = pct_reads_promoter
        v[13] = pct_reads_peaks
        v[14] = pct_mito_reads
        # convert to string for writing out to file
        v = [str(x) for x in v]
        metadata_dict[k] = v
        
    return(metadata_dict)

def write_barcode_metadata(metadata_dict, output_file):
    # write out all fields except for TSS enrichment (must be calculated separately)
    fields = ["barcode","fragments_promoter","reads_tss","reads_promoter","reads_tss_total",
              "reads_unique","reads_duplicate","pct_duplicates","reads_peaks","fragment_peaks",
              "raw_reads_nonmito","raw_reads_mito","pct_reads_promoter","pct_reads_peaks","pct_mito_reads"]

    with open(output_file, "w") as f:
        # write header
        f.write("\t".join(fields) + "\n")
        # write rows
        for k,v in metadata_dict.items():
            f.write(k + "\t" + "\t".join(v[:4]+v[5:]) + "\n")
               
def main():
    args = parse_arguments()
    output_file = getattr(args, "output_file")
    metadata_files = getattr(args, "metadata_files")
    
    metadata_dict = merge_barcode_metadata(metadata_files)
    
    write_barcode_metadata(metadata_dict, output_file)

if __name__ == "__main__":
    main()
