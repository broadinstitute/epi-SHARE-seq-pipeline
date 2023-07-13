#!/usr/bin/env python3
# coding: utf-8
"""
Merge RNA count matrices
"""

import argparse
import gzip
import h5py
import numpy as np
import os.path
import tarfile
from collections import defaultdict
from scipy.sparse import csc_matrix

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate a merged h5 count matrix of genes x barcodes")
    parser.add_argument("prefix", help="Prefix for naming output files")
    parser.add_argument("tar_files", nargs="*", help="File names for tar archives, one per matrix to be merged. Each must contain matrix.mtx, features.tsv, and barcodes.tsv files")
    parser.add_argument("--pkrs", nargs="*", help="PKR names, one per matrix to be merged")
    parser.add_argument("--ensembl", help="Flag for outputting gene names in ENSEMBL form, rather than gene id", action="store_true")
    return parser.parse_args()

def get_split_lines(file_name, delimiter, skip=0):
    """
    Read file contents and yield generator with line entries
    """
    opener = gzip.open if file_name.endswith('.gz') else open
    with opener(file_name, "rt") as f:
        for i in range(skip):
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)
            
def merge_sum_dicts(dicts):
    """
    Merge dictionaries, summing values of repeated keys
    """
    merged = defaultdict(int)
    for d in dicts:
        for k, v in d.items():
            merged[k] += v
    return merged
            
def get_merged_data(tar_files, pkrs, ensembl):
    """
    Takes in paths to tar files (one per count matrix to be merged)
    containing barcodes.tsv, features.tsv, and matrix.mtx output files from STARsolo. 
    Merges matrices by translating individual mtx file indices to a master set of indices. 
    Outputs merged CSC matrix, barcode list, and gene list. 
    If the same barcode(+PKR) is present in multiple tar files, its counts will be summed by default.
    However, if concat==True, each occurrence will be preserved in its own column.
    Genes will be outputted as ENSEMBL IDs if ensembl==True, and as gene symbols otherwise.
    """
    # if PKRs supplied, associate PKRs with appropriate tars
    tar_pkr = dict(zip(tar_files, pkrs)) if pkrs else {}
    gene_set = set()
    barcode_list = []
    count_mappings = []
    barcode_mappings = []
    gene_mappings = []
    ensembl_to_gene = {}
    
    # read in contents of tar files
    for tar_file in tar_files:
        basename = os.path.basename(tar_file)
        tar = tarfile.open(tar_file, mode="r")
        tar.extractall()
        
        # read barcodes file
        barcodes = get_split_lines("barcodes.tsv.gz", delimiter="\t")
        # get mapping of barcode to column index; {barcode:col_idx}
        # append PKR with underscore if supplied
        if tar_pkr:
            barcode_mapping = {(line[0]+"_"+tar_pkr[tar_file]+"_"+basename):idx for idx, line in enumerate(barcodes)}
        else:
            barcode_mapping = {(line[0]+"_"+basename):idx for idx, line in enumerate(barcodes)}
        barcode_list.append(barcode_mapping.keys())
        barcode_mappings.append(barcode_mapping)
        
        # read features file
        features = get_split_lines("features.tsv.gz", delimiter="\t")              
        # get mapping of gene to row index
        if ensembl:
            gene_mapping = {line[0]:idx for idx, line in enumerate(features)}
        else:
            ensembl_dict = {line[0]:line[1] for line in features}
            ensembl_to_gene.update(ensembl_dict)
            gene_mapping = {key:idx for idx, key in enumerate(ensembl_dict.keys())}
            
        gene_mapping = {line[0]:idx for idx, line in enumerate(features)}
        gene_set.update(gene_mapping.keys())
        gene_mappings.append(gene_mapping)
        
        # read matrix file
        matrix = get_split_lines("matrix.mtx.gz", delimiter=" ", skip=3)
        # get mapping of coordinates to count and convert to zero-based indexing;
        # {(row_idx,col_idx):count}
        count_mapping = {(int(line[0])-1, int(line[1])-1) : int(line[2]) for line in matrix}
        count_mappings.append(count_mapping)
        
        tar.close()
    
    barcode_list = list(set(barcode_list))
    gene_list = list(gene_set)
    n_row = len(gene_list)
    n_col = len(barcode_list)
    
    # assign column indices for master list of barcodes    
    merged_barcode_mapping = {barcode:idx for idx, barcode in enumerate(barcode_list)}
    # assign row indices for master list of genes
    merged_gene_mapping = {gene:idx for idx, gene in enumerate(gene_list)}
    
    # iterate through each matrix's mappings and translate them to the merged mappings; 
    # for each gene/barcode, get its index in the merged mapping to create the dict {individual_idx:merged_idx}.
    # then re-map counts into the merged matrix; {(merged_row_idx, merged_col_idx):count}. 
    # make csc matrix.
    count_matrices = list()
    for i in range(len(tar_files)):
        barcode_mappings[i] = {idx:merged_barcode_mapping[barcode] for barcode, idx in barcode_mappings[i].items()}
        gene_mappings[i] = {idx:merged_gene_mapping[gene] for gene, idx in gene_mappings[i].items()}
        count_mappings[i] = {(gene_mappings[i][coords[0]],barcode_mappings[i][coords[1]]):count for coords, count in count_mappings[i].items()}
        count_matrix = csc_matrix((list(count_mappings[i].values()), zip(*count_mappings[i].keys())), shape=(n_row,n_col))
        count_matrices.append(count_matrix)
    
    # sum matrices to get merged matrix
    merged_matrix = csc_matrix((n_row,n_col))
    for matrix in count_matrices:
        merged_matrix += matrix
        
    ## TO DO: collapse rows if using gene name
    
    return(merged_matrix, barcode_list, gene_list)

def write_barcode_metadata(prefix, count_matrix, barcode_list):
    umis = count_matrix.sum(axis=0)
    genes = np.diff(count_matrix.indptr)
    fields = ["barcode", "umis", "genes"]
    
    with open(prefix + "_rna_barcode_metadata.tsv", "w") as f:
        # write header
        f.write("\t".join(fields) + "\n")
        # write rows
        for triple in zip(barcode_list, umis, genes):
            f.write("\t".join(triple[:]) + "\n")
        
def write_h5(prefix, count_matrix, barcode_list, gene_list):
    h5_file = h5py.File(prefix + ".h5", "w")

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
    # get arguments
    args = parse_arguments()
    prefix = getattr(args, "prefix")
    tar_files = getattr(args, "tar_files")
    pkrs = getattr(args, "pkrs")
    ensembl = getattr(args, "ensembl")
    
    # get merged count matrix, barcode list, and gene list from input tars
    count_matrix, barcode_list, gene_list = get_merged_data(tar_files, pkrs, ensembl)
    
    # write barcode metadata file
    write_barcode_metadata(prefix, count_matrix, barcode_list)
    
    # write merged data to new h5 file
    write_h5(prefix, count_matrix, barcode_list, gene_list)

if __name__ == "__main__":
    main()

