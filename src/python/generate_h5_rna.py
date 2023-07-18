#!/usr/bin/env python3
# coding=utf8

"""
This script takes in the STARsolo barcodes tsv file, features tsv file,
and raw count matrix mtx file, and generates an h5 file containing the
genes x barcodes count matrix.
"""

import argparse
from collections import defaultdict
import gzip
import h5py
import logging
from scipy.sparse import csc_matrix

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate an h5 count matrix of genes x barcodes")
    parser.add_argument("matrix_file", help="Filename for STARsolo raw matrix mtx file")
    parser.add_argument("features_file", help="Filename for STARsolo features tsv file")
    parser.add_argument("barcodes_file", help="Filename for STARsolo barcodes tsv file")
    parser.add_argument("output_file", help="Filename for output h5 file")
    parser.add_argument("pkr", help="Experiment prefix", nargs = '?')
    parser.add_argument("--ensembl", help="Flag for outputting genes using ENSEMBL ID, rather than gene name", action="store_true")

    return parser.parse_args()

def get_split_lines(file_name, delimiter, skip=0):
    """Read file contents and yield generator with line entries"""
    opener = gzip.open if file_name.endswith('.gz') else open

    with opener(file_name, "rt") as f:
        for i in range(skip):
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)

def rename_duplicates(duplicate_list):
    """Rename duplicate entries as entry, entry.1, entry.2, etc."""
    seen = defaultdict(int)
    renamed_list = []
    
    for entry in duplicate_list:
        renamed_list.append(f"{entry}.{seen[entry]}" if entry in seen else entry)
        seen[entry] += 1
        
    return renamed_list

def build_count_matrix(matrix_file):
    """Convert contents of mtx file to csc matrix"""
    matrix = get_split_lines(matrix_file, delimiter=" ", skip=2)
    dimensions = next(matrix)
    n_rows = int(dimensions[0])
    n_cols = int(dimensions[1])
    
    count_mapping = {}
    
    for entry in matrix:
        # subtract 1 from indices to convert to zero-based indexing
        row_ind = int(entry[0])-1
        col_ind = int(entry[1])-1
        count = int(entry[2])
        count_mapping[(row_ind,col_ind)] = count
                
    count_matrix = csc_matrix((list(count_mapping.values()), zip(*count_mapping.keys())), shape=(n_rows,n_cols))

    return count_matrix

def write_h5(output_file, count_matrix, barcode_list, gene_list):
    h5_file = h5py.File(output_file, "w")

    # create datasets expected for Seurat import
    g = h5_file.create_group("group")
    g.create_dataset("barcodes", data=barcode_list)
    g.create_dataset("data", data=count_matrix.data)
    g.create_dataset("gene_names", data=gene_list)
    g.create_dataset("genes", data=gene_list)
    g.create_dataset("indices", data=count_matrix.indices)
    g.create_dataset("indptr", data=count_matrix.indptr)
    g.create_dataset("shape", data=count_matrix.shape)

    h5_file.close()

def main():
    # create log file
    logging.basicConfig(filename="generate_h5_rna.log", level=logging.INFO)

    # get arguments
    args = parse_arguments()
    matrix_file = getattr(args, "matrix_file")
    features_file = getattr(args, "features_file")
    barcodes_file = getattr(args, "barcodes_file")
    pkr = getattr(args, "pkr", None)
    output_file = getattr(args, "output_file")
    ensembl = getattr(args, "ensembl")

    # read input files
    logging.info("Reading input files\n")
    
    # get genes from features file
    features = get_split_lines(features_file, delimiter="\t")
    if ensembl:
        gene_list = [line[0] for line in features]
    else:
        gene_list_duplicated = [line[1] for line in features]
        # append .1, .2, etc. for duplicated genes
        gene_list = rename_duplicates(gene_list_duplicated)

    # get barcodes from barcodes file, reformat as R1R2R3_PKR
    barcodes = get_split_lines(barcodes_file, delimiter="\t")
    barcode_list = [line[0] for line in barcodes]
    if pkr is None:
        formatted_barcode_list = barcode_list
    else:
        formatted_barcode_list = [barcode + "_" + pkr for barcode in barcode_list]

    # generate count matrix
    logging.info("Generating count matrix\n")
    count_matrix = build_count_matrix(matrix_file)

    # write h5 file
    logging.info(f"Writing to {output_file}.h5\n")
    write_h5(output_file, count_matrix, formatted_barcode_list, gene_list)
    logging.info("Finished writing h5 file\n")

if __name__ == "__main__":
    main()
