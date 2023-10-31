"""
This script takes the downloaded h5ad file as input, and
outputs an h5ad file which includes the metacell profiles.

The metacells are created using k-nearest-neighbor apparoach
based on stratified-sampling results.
"""
import os
import argparse
import logging
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import scipy as sp
from sklearn.neighbors import NearestNeighbors


def parse_arguments():
    parser = argparse.ArgumentParser(description="Create metacells for reference dataset")
    parser.add_argument("--input", type=str, required=True,
                        help="Filename dataset id to download.", default='reference')
    parser.add_argument("--reference_label", type=str, default='cell_type',
                        help="Category used to down sample the cells. Usually set to cell_type")
    parser.add_argument("--frac", type=float, default=1.0,
                        help="Fraction of cells per category for sub-sampling. Must be between 0 and 1.")
    parser.add_argument("--n_neighbors", type=int, default=10,
                        help="Number of neighbors for aggregation. Default: 10")
        
    return parser.parse_args()


def main():
    # get arguments
    args = parse_arguments()
    
    # create log file
    logging.basicConfig(filename="create_metacells.log", level=logging.INFO)
    
    adata = sc.read_h5ad(f"{args.input}.h5ad")
    
    adata.layers['counts'] = adata.X.copy()
    
    # if X_pca does not exist, running the preprocessing steps to 
    # generate X_pca
    if not 'X_pca' in adata.obsm:
        logging.info("Data processing!")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        sc.tl.pca(adata)
    
    
    # stratified sub-sampling
    obs = adata.obs.groupby(args.reference_label, group_keys=False).apply(lambda x: x.sample(frac=args.frac))
    adata_sub = adata[obs.index].copy()
    
    nbrs = NearestNeighbors(n_neighbors=args.n_neighbors).fit(adata.obsm['X_pca'])
    _, indices = nbrs.kneighbors(adata_sub.obsm['X_pca'])
    
    
    # create metacell profiles
    meta_counts = np.empty(shape=(adata_sub.n_obs, adata_sub.n_vars), dtype=np.float32)
    for i in range(indices.shape[0]):
        meta_counts[i] = np.sum(adata.layers['counts'][indices[i]], axis=0)

    meta_counts = sp.sparse.csr_matrix(meta_counts)
    adata_meta = ad.AnnData(X=meta_counts, 
                            obs=adata_sub.obs,
                           var=adata_sub.var)
    adata_meta.raw = adata_meta.copy()

    adata_meta.write_h5ad(f"{args.input}_meta.h5ad")
    
    logging.info("All done!")


if __name__ == '__main__':
    main()
    
