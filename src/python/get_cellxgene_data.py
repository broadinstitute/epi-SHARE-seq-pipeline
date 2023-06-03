"""
This script takes dataset id as input and download an h5ad file 
from cellxgene server using cellxgene_census API.
"""

import argparse
import logging
import cellxgene_census
import scanpy as sc

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download data from cellxgene server")
    parser.add_argument("--id", type=str, required=True,
                        help="Cellxgene dataset id to download.")
    parser.add_argument("--out", type=str, required=True,
                        help="Output filename", default="reference")
    
    return parser.parse_args()


if __name__ == '__main__':
    # create log file
    logging.basicConfig(filename="get_cellxgene_data.log", level=logging.INFO)

    # get arguments
    args = parse_arguments()
    
    logging.info("Downloading data\n")
    cellxgene_census.download_source_h5ad(
        dataset_id=args.id, 
        to_path=f"{args.out}.h5ad")
    
    adata = sc.read_h5ad(f"{args.out}.h5ad")
    
    logging.info("All done!")