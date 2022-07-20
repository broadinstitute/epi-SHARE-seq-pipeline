#!/usr/bin/env python3

"""
This script takes in a file containing genes, barcodes, and unique read counts,
and generates a genes x barcodes count matrix.
The matrix is outputted as an h5 file.
"""

import argparse
import h5py
import numpy as np
from scipy.sparse import csc_matrix

parser = argparse.ArgumentParser(description="Generate an h5 count matrix of genes x barcodes")
parser.add_argument("-i", "--input_file", help="Input file containing barcodes, genes, unique counts")
parser.add_argument("-o", "--output_prefix", help="Prefix for naming h5 output file")

args = parser.parse_args()

if getattr(args, "input_file") is None:
    print("ERROR: Input file not provided\n")
    parser.parse_args(["-h"])

input_file = getattr(args, "input_file")
output_prefix = getattr(args, "output_prefix")

# use generator to split lines rather than reading into memory
def get_split_lines(file_name):
    with open(file_name, "r") as f:
        for line in f:
            yield line.rstrip().split(sep="\t")

# read in file and build count matrix            
def build_count_matrix(file_name):
    unique_barcodes = {line[0] for line in get_split_lines(file_name)}
    # assign each barcode a column number
    barcode_index_mappings = {barcode:idx for idx, barcode in enumerate(unique_barcodes)}
    # get vector of column indices for input data
    barcode_indices = [barcode_index_mappings[line[0]] for line in get_split_lines(file_name)]

    unique_genes = {line[1] for line in get_split_lines(file_name)}
    # assign each gene a row number 
    gene_index_mappings = {gene:idx for idx, gene in enumerate(unique_genes)}
    # get vector of row indices for input data
    gene_indices = [gene_index_mappings[line[1]] for line in get_split_lines(file_name)]
    
    counts = [int(line[2]) for line in get_split_lines(file_name)]
    
    n_rows = len(unique_genes)
    n_cols = len(unique_barcodes)
    
    # need to use csc for indices to be compatible with downstream Seurat import
    # build matrix by providing row and column indices of each nonzero entry
    count_matrix = csc_matrix((counts, (gene_indices,barcode_indices)), shape=(n_rows,n_cols), dtype=np.int32)
    
    return count_matrix, unique_barcodes, unique_genes

count_matrix, unique_barcodes, unique_genes = build_count_matrix(input_file)

print(f"Number of genes x barcodes = {len(unique_barcodes)} x {len(unique_genes)}\n")

# create h5 file
if getattr(args, "output_prefix") is None:
    print("h5 output prefix not provided; writing to out.gene.bc.matrices.h5\n")
    h5_file = h5py.File("out.gene.bc.matrices.h5", "w")

else: 
    print(f"Writing to {output_prefix}.h5\n")
    h5_file = h5py.File((output_prefix+".h5"), "w")

# create datasets expected for Seurat import    
g = h5_file.create_group("group")
g.create_dataset("barcodes", data=list(unique_barcodes))
g.create_dataset("data", data=count_matrix.data)
g.create_dataset("gene_names", data=list(unique_genes))
g.create_dataset("genes", data=list(unique_genes))
g.create_dataset("indices", data=count_matrix.indices)
g.create_dataset("indptr", data=count_matrix.indptr)
g.create_dataset("shape", data=count_matrix.shape)
h5_file.close()

print("Finished writing h5 file\n")
