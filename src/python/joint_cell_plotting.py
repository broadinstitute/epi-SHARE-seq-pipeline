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
    parser.add_argument("rna_metrics_file", help="Filename for RNA metrics tsv file")
    parser.add_argument("atac_metrics_file", help="Filename for ATAC metrics tsv file")
    parser.add_argument("remove_low_yielding_cells", type=int, help="Minimum number of UMIs/fragments required for a cell to be plotted")
    parser.add_argument("min_umis", type=int, help="Cutoff for minimum number of UMIs")
    parser.add_argument("min_genes", type=int, help="Cutoff for minimum number of genes")
    parser.add_argument("min_tss", type=int, help="Cutoff for minimum TSS score")
    parser.add_argument("min_frags", type=int, help="Cutoff for minimum number of ATAC fragments")
    parser.add_argument("plot_file", help="Filename for plot png file")
    parser.add_argument("qc_summary_data", help="text file that numbers for the final pipeline output are written to")
    parser.add_argument("barcode_metadata_file", help="Filename for barcode metadata csv file")
    parser.add_argument("pkr", help="PKR name", nargs='?', default="")

    return parser.parse_args()

def get_split_lines(file_name, delimiter, skip_header):
    with open(file_name, "r") as f:
        if skip_header:
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)

def merge_dicts(dict_1, dict_2):
    """Merge dictionaries by key; combine values into quadruple, fill with 0s if key not in both dicts"""
    keys = set(dict_1.keys() | dict_2.keys())
    merged = {k: (dict_1.get(k, (0,0)) + dict_2.get(k, (0,0))) for k in keys}

    return(merged)

def get_metrics(rna_metrics_file, atac_metrics_file, remove_low_yielding_cells):
    """Read files and aggregate metrics into Pandas dataframe"""
    rna_metrics_contents = get_split_lines(rna_metrics_file, delimiter="\t", skip_header=True)
    umis = []
    genes = []
    rna_barcodes = []
    # remove cells that have fewer than 10 UMIs
    for line in rna_metrics_contents:
        if int(line[3]) >= remove_low_yielding_cells:
            umis.append(int(line[3]))
            genes.append(int(line[4]))
            rna_barcodes.append(line[0])
    rna_metrics = dict(zip(rna_barcodes, zip(umis, genes)))

    atac_metrics_contents = get_split_lines(atac_metrics_file, delimiter="\t", skip_header=True)
    tss = []
    frags = []
    atac_barcodes = []
    # remove cells that have fewer than 10 fragments
    for line in atac_metrics_contents:
        if int(line[6])/2 >= remove_low_yielding_cells:
            tss.append(float(line[4]))
            frags.append(int(line[6])/2)
            atac_barcodes.append(line[0])
    atac_metrics = dict(zip(atac_barcodes, zip(tss, frags)))

    # merge metrics by barcodes
    metrics = merge_dicts(rna_metrics, atac_metrics)
    df = pd.DataFrame.from_dict(metrics, orient="index", columns=["umis","genes","tss","frags"])

    return(df)

def qc_cells(df, min_umis, min_genes, min_tss, min_frags):
    pass_umis = df["umis"] >= min_umis
    pass_genes = df["genes"] >= min_genes
    pass_tss = df["tss"] >= min_tss
    pass_frags = df["frags"] >= min_frags

    # add df column with QC outcome
    qc_conditions  = [(pass_umis & pass_genes & pass_tss & pass_frags),
                      (pass_umis & pass_genes),
                      (pass_tss & pass_frags),
                      (~(pass_umis & pass_genes) & (~(pass_tss & pass_frags)))]
    qc_choices = ["both", "RNA only", "ATAC only", "neither"]
    df["QC"] = np.select(qc_conditions, qc_choices)

    # get counts of each outcome type (used in plot legend)
    outcome_counts = df["QC"].value_counts()

    df["QC_count"] = [f"{outcome} ({outcome_counts[outcome]})" for outcome in df["QC"]]

    return(df)

def round_to_power_10(x):
    return(10**np.ceil(np.log10(x)))

def label_func(breaks):
    return [int(x) for x in breaks]

def plot_cells(df, pkr, min_umis, min_genes, min_tss, min_frags, plot_file):
    # get max x and y coords to set plot limits
    max_x = max(df["frags"])
    max_y = max(df["umis"])
    xy_lim = round_to_power_10(max(max_x, max_y))

    plot = (ggplot(df, aes("frags", "umis", color="QC_count"))
             + geom_point(size=0.5)
             + labs(title = f"Joint Cell Calling ({pkr})",
                    caption = f"ATAC cutoffs: TSS ≥ {min_tss}, frags ≥ {min_frags}. RNA cutoffs: UMIs ≥ {min_umis}, genes ≥ {min_genes}",
                    x = "ATAC Unique Fragments per Barcode",
                    y = "RNA UMIs per Barcode",
                    color = "QC")
             + theme_light()
             + theme(figure_size = (8,6),
                     title = element_text(size=12),
                     axis_title = element_text(size=10),
                     axis_text = element_text(size=8),
                     legend_box_margin = 0,
                     legend_title = element_text(size=8),
                     legend_text = element_text(size=6),
                     legend_key = element_blank(),
                     plot_caption=element_text(size=8, ha="center", margin={"r": 3.2, "t": -0.2, "units": "in"}),
                     panel_grid_minor = element_blank())
             + scale_x_log10(limits=(10,xy_lim), labels=label_func)
             + scale_y_log10(limits=(10,xy_lim), labels=label_func)
             )

    plot.save(filename=plot_file, dpi=1000)

def write_top_level_txt(input_file, output_file, min_tss, min_frags, min_umis, min_genes): 
    """write the counts of atac, rna, both, and neither to a text file so that they can be 
       outputs of the pipeline"""


    #open file that has the values for the reads of atac, rna, both, and neither
    data = pd.read_csv(input_file)
    
    #get counts and cast to strings so they can be written to the file
    neither_count = get_count_of_type(data, 'neither')
    neither_count_str = str(neither_count) + '\n'
    both_count = get_count_of_type(data, 'both')
    both_count_str = str(both_count) + '\n'
    rna_count = get_count_of_type(data, 'RNA only')
    rna_count_str = str(rna_count) + '\n'
    atac_count = get_count_of_type(data, 'ATAC only')
    atac_count_str = str(atac_count) + '\n'
    
    #cast passed in numbers for thresholds to strings so they can be writen 
    #to the file
    min_tss_str = str(min_tss) + '\n'
    min_frags_str = str(min_frags) + '\n'
    min_umis_str = str(min_umis) + '\n'
    min_genes_str = str(min_genes) 
    
    #open the specified output file and write all the values
    output_file = open(output_file, 'w')
    output_file.write(neither_count_str)
    output_file.write(both_count_str)
    output_file.write(rna_count_str)
    output_file.write(atac_count_str)
    output_file.write(min_tss_str)
    output_file.write(min_frags_str)
    output_file.write(min_umis_str)
    output_file.write(min_genes_str)


def get_count_of_type(dataframe, classifier):
    """extract the count information from the QC column for a given category (both, neither, rna only, atac only)"""
    
    data_one_kind = dataframe.loc[dataframe['QC'] == classifier]
    #return with the value 0 if QC column is never equal to the specified
    #classifier
    if data_one_kind.empty:
        return 0
    
    #extract number contents of QC_count field of the first entry (QC count 
    # should be the same for any of the entries in the dataframe)
    field = data_one_kind.at[data_one_kind.index[1], 'QC_count']
    count = ""
    # extract the digits from other characters in the string
    for char in field: 
        if char.isdigit():
            count = count + char
    return int(count)


def main():
    # create log file
    logging.basicConfig(filename="joint_cell_plotting.log", level=logging.INFO)

    # get arguments
    args = parse_arguments()
    pkr = getattr(args, "pkr")
    rna_metrics_file = getattr(args, "rna_metrics_file")
    atac_metrics_file = getattr(args, "atac_metrics_file")
    remove_low_yielding_cells = getattr(args, "remove_low_yielding_cells")
    barcode_metadata_file = getattr(args, "barcode_metadata_file")
    min_umis = getattr(args, "min_umis")
    min_genes = getattr(args, "min_genes")
    min_tss = getattr(args, "min_tss")
    min_frags = getattr(args, "min_frags")
    plot_file = getattr(args, "plot_file")
    qc_summary_data = getattr(args, "qc_summary_data")

    # read rna and atac files, get cell metrics
    logging.info("Getting metrics\n")
    metrics_df = get_metrics(rna_metrics_file, atac_metrics_file, remove_low_yielding_cells)

    # QC cells based on inputted cutoffs
    logging.info("QCing cells\n")
    metrics_df = qc_cells(metrics_df, min_umis, min_genes, min_tss, min_frags)

    # generate plot
    logging.info("Generating joint cell calling plot\n")
    plot_cells(metrics_df, pkr, min_umis, min_genes, min_tss, min_frags, plot_file)

    # write the stats for the top level into a csv, containd joint qc numbers
    write_top_level_txt(barcode_metadata_file, qc_summary_data, min_tss, min_frags, min_umis, min_genes)

    
    # save dataframe
    logging.info("Saving dataframe as csv\n")
    metrics_df.to_csv(barcode_metadata_file)
    logging.info("All done!")


if __name__ == "__main__":
    main()

