#!/usr/bin/env python3
# coding=utf8

"""
This script takes in the STARsolo barcodes tsv file, features tsv file,
and raw count matrix mtx file, and generates an h5 file containing the
genes x barcodes count matrix.
"""

import argparse
import gzip
import h5py
import logging
from scipy.sparse import csc_matrix

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate an h5 count matrix of genes x barcodes")
    parser.add_argument("matrix_file", help="Filename for STARsolo raw matrix mtx file")
    parser.add_argument("features_file", help="Filename for STARsolo features tsv file")
    parser.add_argument("barcodes_file", help="Filename for STARsolo barcodes tsv file")
    parser.add_argument("pkr", help="Experiment prefix")
    parser.add_argument("output_file", help="Filename for output h5 file")

    return parser.parse_args()

def get_split_lines(file_name, delimiter, skip=0):
    """Read file contents and yield generator with line entries"""
    opener = gzip.open if file_name.endswith('.gz') else open

    with opener(file_name, "rt") as f:
        for i in range(skip):
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)

def build_count_matrix(matrix):
    """Convert contents of mtx file to csc matrix"""
    # first line of matrix contains dimensions
    dimensions = next(matrix)
    n_rows = int(dimensions[0])
    n_cols = int(dimensions[1])

    gene_indices = []
    barcode_indices = []
    counts = []

    for line in matrix:
        # subtract 1 from indices to convert to zero-based indexing
        gene_indices.append(int(line[0])-1)
        barcode_indices.append(int(line[1])-1)
        counts.append(int(line[2]))

    count_matrix = csc_matrix((counts, (gene_indices,barcode_indices)), shape=(n_rows,n_cols))

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
    pkr = getattr(args, "pkr")
    output_file = getattr(args, "output_file")

    # read input files
    logging.info("Reading input files\n")
    # skip first two lines of matrix file (header)
    matrix = get_split_lines(matrix_file, delimiter=" ", skip=2)
    # get genes from features file
    features = get_split_lines(features_file, delimiter="\t")
    gene_list = [line[1] for line in features]
    # get barcodes from barcodes file, reformat as R1,R2,R3,PKR
    barcodes = get_split_lines(barcodes_file, delimiter="\t")
    barcode_list = [line[0] for line in barcodes]
    formatted_barcode_list = [barcode[:8] + "," + barcode[8:16] + "," + barcode[16:] + "," + pkr for barcode in barcode_list]

    # generate count matrix
    logging.info("Generating count matrix\n")
    count_matrix = build_count_matrix(matrix)

    # write h5 file
    logging.info(f"Writing to {output_file}.h5\n")
    write_h5(output_file, count_matrix, formatted_barcode_list, gene_list)
    logging.info("Finished writing h5 file\n")

if __name__ == "__main__":
    main()
