#!/usr/bin/env python3

"""
This script takes in an RNA bed file containing read counts and generates a count matrix of genes x cells.
The matrix is outputted both as a csv for downstream QC/plotting and as an h5ad file.
"""

import anndata as ad
import argparse
import pandas as pd
from scipy.sparse import csr_matrix

parser = argparse.ArgumentParser(description="Generate an h5 count matrix of genes x cells")
parser.add_argument("-b", help="bed file name")

args = parser.parse_args()
bedfile = args.b

print("Loading linear gene table...")
linear = pd.read_csv(bedfile, delimiter="\t", names=["barcode","gene","unique_counts","dup_counts"])
print("Finished loading file")

cells = set(linear["barcode"])
genes = set(linear["gene"])

print(f"No. genes x cells = {len(genes)} x {len(cells)}")

# convert bed to genes x cell dataframe
count_df = linear.pivot(index="gene", columns="barcode", values="unique_counts")
count_df = count_df.fillna(0)

# convert count df to compressed sparse row matrix
count_mat = csr_matrix(count_df)

# convert sparse matrix to AnnData object
count_ad = ad.AnnData(count_mat)

# add gene names as observations (rownames), add cell barcodes as variables (colnames)
count_ad.obs_names = count_df.index
count_ad.var_names = count_df.columns

print("Writing to out.gene.bc.matrices.csv")
count_df.to_csv("out.gene.bc.matrices.csv", index=False)
print("Finished writing csv file")

print("Writing to out.gene.bc.matrices.h5ad")
count_ad.write(filename="out.gene.bc.matrices.h5ad")
print("Finished writing h5ad file")
