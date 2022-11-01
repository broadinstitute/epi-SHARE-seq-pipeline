#!/usr/bin/env python3

"""
This script QCs barcodes via ATAC frags & TSS and RNA UMIs & genes,
and plots all barcodes colored by joint QC status. It also generates the 
same plot with transparency added to show density. 
"""

import argparse
import logging
import numpy as np
import pandas as pd
from plotnine import *

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot barcodes by RNA and ATAC QC status")
    parser.add_argument("pkr", help="PKR name")
    parser.add_argument("rna_metrics_file", help="Filename for RNA metrics tsv matrix")
    parser.add_argument("atac_metrics_file", help="Filename for ATAC metrics tsv file")
    parser.add_argument("min_UMIs", type=int, help="Cutoff for minimum number of UMIs")
    parser.add_argument("min_genes", type=int, help="Cutoff for minimum number of genes")
    parser.add_argument("min_TSS", type=int, help="Cutoff for minimum TSS score")
    parser.add_argument("min_frags", type=int, help="Cutoff for minimum number of ATAC fragments")
    parser.add_argument("barcode_metadata_file", help="Filename for barcode metadata csv file")
    
    return parser.parse_args()

def get_split_lines(file_name, delimiter):
    with open(file_name, "r") as f:
        for line in f:
            yield line.rstrip().split(sep=delimiter)

def merge_dicts(dict_1, dict_2):
    """Merge dictionaries by key; combine values into quadruple, fill with 0s if key not in both dicts"""
    keys = set(dict_1.keys() | dict_2.keys())
    merged = {k: (dict_1.get(k, (0,0)) + dict_2.get(k, (0,0))) for k in keys}
    
    return(merged)

def get_metrics(rna_metrics_file, atac_metrics_file):
    """Read files and aggregate metrics into Pandas dataframe"""
    rna_metrics_contents = get_split_lines(rna_metrics_file, delimiter="\t")
    UMIs = []
    genes = []
    rna_barcodes = []
    next(rna_metrics_contents) # skip header
    for line in rna_metrics_contents:
        if int(line[1]) >= 10: # impose cutoff of 10 UMIs; RNA metrics file contains all values
            UMIs.append(int(line[1]))
            genes.append(int(line[2]))
            rna_barcodes.append(line[5])
    rna_metrics = dict(zip(rna_barcodes, zip(UMIs, genes)))
    
    atac_metrics_contents = get_split_lines(atac_metrics_file, delimiter="\t")
    TSS = []
    frags = []
    atac_barcodes = []
    next(atac_metrics_contents) # skip header 
    for line in atac_metrics_contents: 
        TSS.append(float(line[1]))
        frags.append(int(line[10]))
        atac_barcodes.append(line[15].split("#")[1]) # get barcode string after library name
    atac_metrics = dict(zip(atac_barcodes, zip(TSS, frags)))
    
    # merge metrics by barcodes
    metrics = merge_dicts(rna_metrics, atac_metrics)
    df = pd.DataFrame.from_dict(metrics, orient="index", columns=["UMIs","genes","TSS","frags"])
    
    return(df)

def qc_cells(df, min_UMIs, min_genes, min_TSS, min_frags):
    u = df["UMIs"] >= min_UMIs
    g = df["genes"] >= min_genes
    t = df["TSS"] >= min_TSS
    f = df["frags"] >= min_frags
    
    # add df column with QC outcome
    qc_conditions  = [(u & g & t & f), (u & g), (t & f)]
    qc_choices = ["both", "RNA only", "ATAC only"]
    df["QC"] = np.select(qc_conditions, qc_choices, default="neither")
    
    # get counts of each outcome type (used in plot legend)
    outcome_counts = df["QC"].value_counts().sort_index() # sort to avoid ordering by count value
    outcome_conditions = [df["QC"]=="ATAC only", df["QC"]=="RNA only", df["QC"]=="both", df["QC"]=="neither"]
    count_choices = [f"{outcome} ({outcome_counts[outcome]})" for outcome in outcome_counts.index]
    df["QC_count"] = np.select(outcome_conditions, count_choices)
    
    return(df)

def round_to_power_10(x):
    return(10**np.ceil(np.log10(x)))

def label_func(breaks):
    return [int(x) for x in breaks]
    
def plot_cells(df, pkr, min_UMIs, min_genes, min_TSS, min_frags):
    # get max x and y coords to set plot limits
    max_x = max(df["frags"])
    max_y = max(df["UMIs"])
    xy_lim = round_to_power_10(max(max_x, max_y))
    
    plot = (ggplot(df, aes("frags", "UMIs", color="QC_count"))
     + geom_point(size=0.5)
     + labs(title = f"Joint Cell Calling ({pkr})",
            caption = f"ATAC cutoffs: TSS ≥ {min_TSS}, frags ≥ {min_frags}. RNA cutoffs: UMIs ≥ {min_UMIs}, genes ≥ {min_genes}",
            x = "ATAC Unique Fragments per Barcode",
            y = "RNA UMIs per Barcode",
            color = "QC")
     + theme_light()
     + theme(figure_size = (8,6),
             legend_box_margin = 0,
             legend_title = element_text(size=10),
             legend_text = element_text(size=8),
             legend_key = element_blank(),
             plot_caption=element_text(size=10, ha="center", margin={"r": 3, "t": -0.1, "units": "in"}),
             panel_grid_minor = element_blank())
     + scale_x_log10(limits=(10,xy_lim), labels=label_func)
     + scale_y_log10(limits=(10,xy_lim), labels=label_func)
     )
    
    plot.save(filename = f"{pkr}_joint_cell_plot.png", dpi=1000)
    
def main():
    # create log file
    logging.basicConfig(filename="joint_cell_plotting.log", level=logging.INFO)
    
    # get arguments
    args = parse_arguments() 
    pkr = getattr(args, "pkr")
    rna_metrics_file = getattr(args, "rna_metrics_file")
    atac_metrics_file = getattr(args, "atac_metrics_file")
    min_UMIs = getattr(args, "min_UMIs")
    min_genes = getattr(args, "min_genes")
    min_TSS = getattr(args, "min_TSS")
    min_frags = getattr(args, "min_frags")
    barcode_metadata_file = getattr(args, "barcode_metadata_file")
    
    # read rna and atac files, get cell metrics
    logging.info("Getting metrics\n")
    metrics_df = get_metrics(rna_metrics_file, atac_metrics_file)
    
    # QC cells based on inputted cutoffs
    logging.info("QCing cells\n")
    metrics_df = qc_cells(metrics_df, min_UMIs, min_genes, min_TSS, min_frags)
    
    # generate plot
    logging.info("Generating joint cell calling plot\n")
    plot_cells(metrics_df, pkr, min_UMIs, min_genes, min_TSS, min_frags)
    
    # save dataframe
    logging.info("Saving dataframe as csv\n")
    metrics_df.to_csv(barcode_metadata_file)
    logging.info("All done!")
    

if __name__ == "__main__":
    main()

