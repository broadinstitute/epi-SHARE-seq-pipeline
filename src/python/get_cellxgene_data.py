#!/usr/bin/env python3
# coding=utf8

"""
This script takes dataset id as input and download an h5ad file 
from cellxgene server using cellxgene_census API.
"""

import argparse
import logging
import cellxgene_census

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download data from cellxgene server")
    parser.add_argument("dataset_id", help="Dataset id to download")
    parser.add_argument("output_filename", help="Output filename")
    
    return parser.parse_args()


if __name__ == '__main__':
    
    # create log file
    logging.basicConfig(filename="get_cellxgene_data.log", level=logging.INFO)

    # get arguments
    args = parse_arguments()
    
    logging.info("Downloading data\n")
    cellxgene_census.download_source_h5ad(
        args.dataset_id, 
        to_path=args.output_filename)
    
    logging.info("All done!")