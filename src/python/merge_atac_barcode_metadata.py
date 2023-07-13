#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge ATAC barcode metadata TSV files
"""

import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge ATAC barcode metadata TSV files")
    parser.add_argument("output_file", help="Filename for output merged barcode metadata TSV file")
    parser.add_argument("metadata_files", nargs="*", help="Filenames for barcode metadata TSV files to be merged")
    return parser.parse_args()
            
def merge_barcode_metadata(metadata_files):
    # only need columns not computed by TSS enrichment script
    columns = ["barcode",
               "reads_unique",
               "reads_duplicate",
               "pct_duplicates",
               "reads_peaks",
               "fragment_peaks",
               "raw_reads_nonmito",
               "raw_reads_mito",
               "pct_reads_promoter",
               "pct_reads_peaks",
               "pct_mito_reads"]
    
    # read metadata into data frames
    metadata = [pd.read_table(metadata_file, header=0, usecols=columns) for metadata_file in metadata_files]
    
    # concatenatate data frames
    merged = pd.concat(metadata)

    # collapse duplicate barcodes by operating on columns accordingly
    agg_dict = {"reads_unique":"sum",
                "reads_duplicate":"sum",
                "pct_duplicates":"mean",
                "reads_peaks":"sum",
                "fragment_peaks":"sum",
                "raw_reads_nonmito":"sum",
                "raw_reads_mito":"sum",
                "pct_reads_promoter":"mean",
                "pct_reads_peaks":"mean",
                "pct_mito_reads":"mean"}    
    merged = merged.groupby("barcode").agg(agg_dict).reset_index()
    
    # sort barcodes to ensure consistency with TSS enrichment tsv
    merged = merged.sort_values("barcode")
        
    return(merged)
    
def main():
    args = parse_arguments()
    output_file = getattr(args, "output_file")
    metadata_files = getattr(args, "metadata_files")
    
    merged = merge_barcode_metadata(metadata_files)
    
    merged.to_csv(output_file, sep="\t", index=False)    

if __name__ == "__main__":
    main()
