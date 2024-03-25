#!/usr/bin/env python3
# coding: utf-8
"""
Take multiple STARsolo mtx tar files and combine to create a merged h5 count matrix.
Also produces a merged version of the STARsolo mtx tar files, as well as a TSV
mapping barcodes to their dataset of origin.
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
    parser.add_argument("merged_h5_file", help="File name for output merged h5 file")
    parser.add_argument("dataset_barcodes_file", help="File name for output dataset barcodes tsv file")
    parser.add_argument("--tar_files", nargs="*", help="File names for tar archives, one per matrix to be merged. Each must contain matrix.mtx, features.tsv, and barcodes.tsv files")
    parser.add_argument("--subpools", nargs="*", help="Cellular sub-pool names, one per matrix to be merged")
    parser.add_argument("--datasets", nargs="*", help="Dataset names, one per matrix to be merged")
    parser.add_argument("--ensembl", help="Flag for outputting gene names in ENSEMBL form, rather than gene id", action="store_true")
    return parser.parse_args()


def get_split_lines(file_name, delimiter, skip=0):
    """
    Read file contents and yield generator with line entries
    """
    opener = gzip.open if file_name.endswith(".gz") else open
    with opener(file_name, "rt") as f:
        for i in range(skip):
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)


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

            
def get_merged_data(tar_files, subpools, datasets, ensembl):
    """
    Takes in paths to tar files (one per count matrix to be merged)
    containing barcodes.tsv, features.tsv, and matrix.mtx output files from STARsolo. 
    Merges matrices by translating individual mtx file indices to a master set of indices. 
    Outputs merged CSC matrix, barcode list, gene list, and barcode-dataset corerspondences. 
    Genes will be outputted as ENSEMBL IDs if ensembl==True, and as gene symbols otherwise.
    """
    mtxs = []
    barcode_mapping_dicts = []
    gene_mapping_dicts = []
    ensembl_to_gene = {}
    barcode_to_dataset = {}
    
    # for each input, read in mtx file, map column indices to barcodes, map row indices to genes
    for i in range(len(tar_files)):
        tar = tarfile.open(tar_files[i], mode="r")
        tar.extract("matrix.mtx.gz")
        tar.extract("barcodes.tsv.gz")
        tar.extract("features.tsv.gz")
        
        # read mtx file
        mtxs.append(get_split_lines("matrix.mtx.gz", delimiter=" ", skip=3))

        # read barcodes file
        barcodes = get_split_lines("barcodes.tsv.gz", delimiter="\t")

        # get mapping of column index to barcode; {col_idx: barcode}
        # append subpool name with underscore if supplied
        if subpools:
            col_to_barcode = {idx: line[0] + "_" + subpools[i] for idx, line in enumerate(barcodes)}
        else:
            basename = os.path.basename(tar_files[i])
            col_to_barcode = {idx: line[0] + "_" + basename for idx, line in enumerate(barcodes)}

        barcode_mapping_dicts.append(col_to_barcode)

        # get mapping of barcodes to dataset name
        for barcode in col_to_barcode.values():
            barcode_to_dataset[barcode] = datasets[i]
        
        # read features file
        features = get_split_lines("features.tsv.gz", delimiter="\t")      

        # make dictionary of {ENSEMBL: gene_id}
        gene_dict = {line[0]: line[1] for line in features}
        ensembl_to_gene.update(gene_dict)

        # get mapping of row index to ENSEMBL id; {row_idx: ENSEMBL}
        row_to_gene = {idx: key for idx, key in enumerate(gene_dict.keys())}
        gene_mapping_dicts.append(row_to_gene)
        
        tar.close()
    
    # get list of all barcodes and all genes
    barcode_list = [k for d in barcode_mapping_dicts for k in d.values()]
    gene_list = [k for d in gene_mapping_dicts for k in d.values()]
    
    # assign column indices to merged barcode list and row indices to merged gene list
    barcode_to_merged_col = {barcode: idx for idx, barcode in enumerate(barcode_list)}
    gene_to_merged_row = {gene: idx for idx, gene in enumerate(gene_list)}
    
    # get mappings of counts to their merged indices
    count_mapping = {}
    for mtx, col_to_barcode, row_to_gene in zip(mtxs, barcode_mapping_dicts, gene_mapping_dicts):
        for entry in mtx:
            row_ind = int(entry[0]) - 1  # indices have 1-based indexing in mtx
            col_ind = int(entry[1]) - 1
            count = int(entry[2])

            merged_row_ind = gene_to_merged_row[row_to_gene[row_ind]]
            merged_col_ind = barcode_to_merged_col[col_to_barcode[col_ind]]

            count_mapping[(merged_row_ind, merged_col_ind)] = count

    # make sparse matrix
    n_row = len(gene_list)
    n_col = len(barcode_list)
    merged_matrix = csc_matrix((list(count_mapping.values()), zip(*count_mapping.keys())), shape=(n_row, n_col))

    # if gene names to be outputted, convert and rename duplicate genes
    if not ensembl:
        gene_list = rename_duplicates([ensembl_to_gene[ensembl_id] for ensembl_id in gene_list])
    
    return(merged_matrix, barcode_list, gene_list, ensembl_to_gene, barcode_to_dataset)

        
def write_h5(count_matrix, barcode_list, gene_list, merged_h5_file):
    h5_file = h5py.File(merged_h5_file, "w")

    g = h5_file.create_group("group")
    g.create_dataset("barcodes", data=barcode_list)
    g.create_dataset("data", data=count_matrix.data)
    g.create_dataset("gene_names", data=gene_list)
    g.create_dataset("genes", data=gene_list)
    g.create_dataset("indices", data=count_matrix.indices)
    g.create_dataset("indptr", data=count_matrix.indptr)
    g.create_dataset("shape", data=count_matrix.shape)

    h5_file.close()   

    
def write_starsolo_outputs(prefix, count_matrix, barcode_list, ensembl_to_gene):
    with gzip.open(prefix + ".features.tsv.gz", "wt") as f:
        for ensembl, gene in ensembl_to_gene.items(): 
            f.write("%s\t%s\n" % (ensembl, gene))
    
    with gzip.open(prefix + ".barcodes.tsv.gz", "wt") as f:
        for barcode in barcode_list:
            f.write("%s\n" % barcode)
            
    row, col = count_matrix.nonzero()
    with gzip.open(prefix + ".matrix.mtx.gz", "wt") as f:
        # write header lines
        f.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        f.write("%s %s %s\n" % (count_matrix.shape[0], count_matrix.shape[1], count_matrix.nnz))
        for triple in zip(row, col, count_matrix.data):
            f.write("%s %s %s\n" % triple)


def write_dataset_barcodes(barcode_to_dataset, dataset_barcodes_file):
    with open(dataset_barcodes_file, "w") as f:
        f.write("barcode\tdataset\n")
        for barcode, dataset in barcode_to_dataset.items():
            f.write(barcode + "\t" + dataset + "\n")	 

            
def main():
    # get arguments
    args = parse_arguments()
    prefix = getattr(args, "prefix")
    merged_h5_file = getattr(args, "merged_h5_file")
    dataset_barcodes_file = getattr(args, "dataset_barcodes_file")
    tar_files = getattr(args, "tar_files")
    subpools = getattr(args, "subpools")
    datasets = getattr(args, "datasets")
    ensembl = getattr(args, "ensembl")
    
    # get merged data from input tars
    count_matrix, barcode_list, gene_list, ensembl_to_gene, barcode_to_dataset = get_merged_data(tar_files, subpools, datasets, ensembl)
    
    # write merged data to h5 file
    write_h5(count_matrix, barcode_list, gene_list, merged_h5_file)
    
    # write merged data to matrix.mtx, features.tsv, and barcodes.tsv files
    write_starsolo_outputs(prefix, count_matrix, barcode_list, ensembl_to_gene)

    # write dataset barcode tsv
    write_dataset_barcodes(barcode_to_dataset, dataset_barcodes_file)

if __name__ == "__main__":
    main()

