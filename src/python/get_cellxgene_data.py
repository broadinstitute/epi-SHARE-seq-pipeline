"""
This script takes dataset id as input and download an h5ad file 
from cellxgene server using cellxgene_census API.
"""

import argparse
import logging
import cellxgene_census
import scanpy as sc
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download data from cellxgene server")
    parser.add_argument("--id", type=str, required=True,
                        help="Cellxgene dataset id to download.")
    parser.add_argument("--out", type=str, required=True,
                        help="Output filename", default="reference")
    parser.add_argument("--reference_label", type=str, default='cell_type',
                        help="Category used to down sample the cells. Usually set to cell_type")
    parser.add_argument("--downsample_frac", type=float, default=1.0,
                        help="Number of cells per category after down sampling.")
        
    return parser.parse_args()


def main():
    # get arguments
    args = parse_arguments()
    
    # create log file
    logging.basicConfig(filename="get_cellxgene_data.log", level=logging.INFO)
    
    logging.info("Downloading data\n")
    
    cellxgene_census.download_source_h5ad(
        dataset_id=args.id, 
        to_path=f"{args.out}.h5ad")
    
    adata = sc.read_h5ad(f"{args.out}.h5ad")
    
    # get counts
    if not adata.raw:
        adata.raw = adata.copy()
    
    # down sample cells stratified by obs_key
    if args.downsample_frac < 1:
        logging.info("Down sampling data\n")
        
        obs = adata.obs.groupby(args.reference_label, group_keys=False).apply(lambda x: x.sample(frac=args.downsample_frac))
        adata = adata[obs.index].copy()

    adata.write_h5ad(f"{args.out}.h5ad")
    
    logging.info("All done!")


if __name__ == '__main__':
    main()
    
