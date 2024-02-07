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
    parser.add_argument("barcode_metadata_file", help="Filename for barcode metadata csv file")
    parser.add_argument("subpool", help="Cellular subpool name", nargs="?", default="")

    return parser.parse_args()


def get_split_lines(file_name, delimiter, skip=0):
    """
    Read file contents and yield generator with line entries
    """
    with open(file_name, "rt") as f:
        for i in range(skip):
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)


def merge_dicts(dict_1, dict_2):
    """Merge dictionaries by key; combine values into quadruple, fill with 0s if key not in both dicts"""
    keys = set(dict_1.keys() | dict_2.keys())
    merged = {k: (dict_1.get(k, (0, 0)) + dict_2.get(k, (0, 0))) for k in keys}

    return merged


def get_metrics(rna_metrics_file, atac_metrics_file, remove_low_yielding_cells):
    """Read files and aggregate metrics into Pandas dataframe"""
    rna_metrics_contents = get_split_lines(rna_metrics_file, delimiter="\t")
    rna_metrics_header = next(rna_metrics_contents)
    rna_barcode_ind = rna_metrics_header.index("barcode")
    umi_ind = rna_metrics_header.index("unique_umi")
    gene_ind = rna_metrics_header.index("genes_final")
    umis = []
    genes = []
    rna_barcodes = []

    for line in rna_metrics_contents:
        if int(line[umi_ind]) >= remove_low_yielding_cells:
            umis.append(int(line[umi_ind]))
            genes.append(int(line[gene_ind]))
            rna_barcodes.append(line[rna_barcode_ind])
    rna_metrics = dict(zip(rna_barcodes, zip(umis, genes)))

    atac_metrics_contents = get_split_lines(atac_metrics_file, delimiter="\t")
    atac_metrics_header = next(atac_metrics_contents)
    atac_barcode_ind = atac_metrics_header.index("barcode")
    reads_ind = atac_metrics_header.index("n_fragment")
    tss_ind = atac_metrics_header.index("tsse")
    tss = []
    frags = []
    atac_barcodes = []
    # remove cells that have fewer than 10 fragments
    for line in atac_metrics_contents:
        if int(line[reads_ind]) >= remove_low_yielding_cells:
            tss.append(float(line[tss_ind]))
            frags.append(int(line[reads_ind]))
            atac_barcodes.append(line[atac_barcode_ind])
    atac_metrics = dict(zip(atac_barcodes, zip(tss, frags)))

    # merge metrics by barcodes
    metrics = merge_dicts(rna_metrics, atac_metrics)
    df = pd.DataFrame.from_dict(metrics, orient="index", columns=["unique_umi", "genes_final", "tss", "frags"])

    return df


def qc_cells(df, min_umis, min_genes, min_tss, min_frags):
    pass_umis = df["unique_umi"] >= min_umis
    pass_genes = df["genes_final"] >= min_genes
    pass_tss = df["tss"] >= min_tss
    pass_frags = df["frags"] >= min_frags

    # add df column with QC outcome
    qc_conditions = [(pass_umis & pass_genes & pass_tss & pass_frags),
                     (pass_umis & pass_genes),
                     (pass_tss & pass_frags),
                     (~(pass_umis & pass_genes) & (~(pass_tss & pass_frags)))]
    qc_choices = ["both", "RNA only", "ATAC only", "neither"]
    df["QC"] = np.select(qc_conditions, qc_choices)

    # get counts of each outcome type (used in plot legend)
    outcome_counts = df["QC"].value_counts()

    df["QC_count"] = [f"{outcome} ({outcome_counts[outcome]})" for outcome in df["QC"]]

    return df


def round_to_power_10(x):
    return 10**np.ceil(np.log10(x))


def label_func(breaks):
    return [int(x) for x in breaks]


def plot_cells(df, subpool, min_umis, min_genes, min_tss, min_frags, plot_file):
    # get max x and y coords to set plot limits
    max_x = max(df["frags"])
    max_y = max(df["unique_umi"])
    xy_lim = round_to_power_10(max(max_x, max_y))

    plot = (ggplot(df, aes("frags", "unique_umi", color="QC_count"))
            + geom_point(size=0.5)
            + labs(title=f"Joint Cell Calling ({subpool})",
                   caption=f"ATAC cutoffs: TSS ≥ {min_tss}, frags ≥ {min_frags}. RNA cutoffs: UMIs ≥ {min_umis}, genes ≥ {min_genes}",
                   x="ATAC Unique Fragments per Barcode",
                   y="RNA UMIs per Barcode",
                   color="QC")
            + theme_light()
            + theme(figure_size=(8, 6),
                    title=element_text(size=12),
                    axis_title=element_text(size=10),
                    axis_text=element_text(size=8),
                    legend_box_margin=0,
                    legend_title=element_text(size=10),
                    legend_text=element_text(size=8),
                    legend_key=element_blank(),
                    plot_caption=element_text(size=8, ha="center", margin={"r": 3.2, "t": 0.1, "units": "in"}),
                    panel_grid_minor=element_blank())
            + scale_x_log10(limits=(10, xy_lim), labels=label_func)
            + scale_y_log10(limits=(10, xy_lim), labels=label_func)
            )

    plot.save(filename=plot_file, dpi=1000)


def main():
    # create log file
    logging.basicConfig(filename="joint_cell_plotting.log", level=logging.INFO)

    # get arguments
    args = parse_arguments()
    subpool = getattr(args, "subpool")
    rna_metrics_file = getattr(args, "rna_metrics_file")
    atac_metrics_file = getattr(args, "atac_metrics_file")
    remove_low_yielding_cells = getattr(args, "remove_low_yielding_cells")
    barcode_metadata_file = getattr(args, "barcode_metadata_file")
    min_umis = getattr(args, "min_umis")
    min_genes = getattr(args, "min_genes")
    min_tss = getattr(args, "min_tss")
    min_frags = getattr(args, "min_frags")
    plot_file = getattr(args, "plot_file")

    # read rna and atac files, get cell metrics
    logging.info("Getting metrics\n")
    metrics_df = get_metrics(rna_metrics_file, atac_metrics_file, remove_low_yielding_cells)

    # QC cells based on inputted cutoffs
    logging.info("QCing cells\n")
    metrics_df = qc_cells(metrics_df, min_umis, min_genes, min_tss, min_frags)

    # generate plot
    logging.info("Generating joint cell calling plot\n")
    plot_cells(metrics_df, subpool, min_umis, min_genes, min_tss, min_frags, plot_file)

    # save dataframe
    logging.info("Saving dataframe as csv\n")
    metrics_df.to_csv(barcode_metadata_file)
    logging.info("All done!")


if __name__ == "__main__":
    main()
