#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge barcode metadata TSV files
"""

import argparse
import pandas as pd
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge barcode metadata TSV files")
    parser.add_argument("modality", choices=["RNA", "ATAC"], help="Sequenced data type")
    parser.add_argument("output_file", help="Filename for output merged barcode metadata TSV file")
    parser.add_argument("metadata_files", nargs="*", help="Filenames for barcode metadata TSV files to be merged")
    parser.add_argument("--concat", help="Flag for appending duplicate barcodes, rather than merging their metadata values", action="store_true")
    return parser.parse_args()
            
def rename_duplicates(duplicate_list):
    """
    Rename duplicate entries as entry, entry.1, entry.2, etc.
    """
    seen = defaultdict(int)
    renamed_list = []
    
    for entry in duplicate_list:
        renamed_list.append(f"{entry}.{seen[entry]}" if entry in seen else entry)
        seen[entry] += 1
        
    return renamed_list
            
def merge_barcode_metadata(metadata_files, modality, concat):
    # read metadata into data frames
    metadata = [pd.read_table(metadata_file, header=0) for metadata_file in metadata_files]
    # concatenatate data frames
    merged = pd.concat(metadata)
    
    if concat:
        # rename duplicate barcodes to make unique
        merged["barcode"] = rename_duplicates(list(merged["barcode"]))
    else:
        # collapse duplicate barcodes by operating on columns accordingly
        if modality == "ATAC":
            agg_dict = {"fragments_promoter":"sum",
                        "reads_tss":"sum",
                        "reads_promoter":"sum",
                        "tss_enrichment":"mean",
                        "reads_tss_total":"sum",
                        "reads_unique":"sum",
                        "reads_duplicate":"sum",
                        "pct_duplicates":"mean",
                        "reads_peaks":"sum",
                        "fragment_peaks":"sum",
                        "raw_reads_nonmito":"sum",
                        "raw_reads_mito":"sum",
                        "pct_reads_promoter":"mean",
                        "pct_reads_peaks":"mean",
                        "pct_mito_reads":"mean"}
        elif modality == "RNA":
            agg_dict = {"total_counts":"sum",
                        "duplicate_counts":"sum",
                        "umis":"sum",
                        "genes":"sum",
                        "percent_mitochondrial":"mean"}
        
        merged = merged.groupby("barcode").agg(agg_dict).reset_index()
        
    return(merged)
    
def main():
    args = parse_arguments()
    modality = getattr(args, "modality")
    output_file = getattr(args, "output_file")
    metadata_files = getattr(args, "metadata_files")
    concat = getattr(args, "concat")
    
    merged = merge_barcode_metadata(metadata_files, modality, concat)
    
    merged.to_csv(output_file, sep="\t", index=False)    

if __name__ == "__main__":
    main()
