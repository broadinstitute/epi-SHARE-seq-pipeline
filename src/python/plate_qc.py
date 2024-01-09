#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:14:50 2023

@author: mknudson
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def parse_arguments():
    parser = argparse.ArgumentParser(description="Make diagram of reads from each round 1 barcode well. Assumes 192 barcodes across two 8x12 plates")
    parser.add_argument("r1_reads_file", help="File containing two tab-separated columns: R1 barcodes, and number of reads per R1 barcode")
    parser.add_argument("r1_barcode_file", help="File containing R1 barcodes, one per line, ordered by well number. Well numbering begins at top left of plate, and proceeds columnwise per plate")
    parser.add_argument("r1_barcode_set_file", help="File containing round 1 barcodes of plate subsets, one subset per line")
    parser.add_argument("prefix", help="Prefix for naming output plot")

    return parser.parse_args()


def get_split_lines(file_name, delimiter, skip=0):
    """Read file contents and yield generator with line entries"""
    with open(file_name, "rt") as f:
        for i in range(skip):
            next(f)
        for line in f:
            yield line.rstrip().split(sep=delimiter)


def get_plate_coords(r1_barcodes, n_rows=8):
    """
    Create dictionary mapping R1 barcodes to their plate coordinates.

    Takes list of R1 barcodes, ordered by well number, and maps each R1 barcode
    to its (col_idx, row_idx) coordinates on an n_rows plate, filled columnwise
    from top left. With two plates, columns will be numbered 0 to 23.
    """
    r1_coords_dict = {}
    for idx, bc in enumerate(r1_barcodes):
        col_idx = idx // n_rows
        row_idx = idx % n_rows
        r1_coords_dict[bc] = (col_idx, row_idx)
    return r1_coords_dict


def get_plate_divisions(r1_barcode_subsets, r1_coords_dict, n_cols=12):
    """
    
    """    
    for subset in r1_barcode_subsets:
        # get column number of each barcode's coordinate tuple
        cols = [r1_coords_dict[bc][0] for bc in subset]
        # get leftmost column of subset
        min_col = min(cols)
        # get rightmost column of subset
        max_col = max(cols)
        
        min_col = min(coords, key=lambda c: c[0])[0]
        max_col = max(coords, key=lambda c: c[0])
        
        max_x = max(tuple_list, key=lambda t: t[0])[0]
max_y = max(tuple_list, key=lambda t: t[1])[1]
        

def plate_plot(r1_reads, r1_barcode_subsets, prefix):
    # https://github.com/karthik/wesanderson
    zissou = ListedColormap(["#3A9AB2", "#449EB5", "#4FA3B8", "#5AA8BB", "#65ADBE",
                             "#6FB2C0", "#76B3BE", "#7DB5BC", "#84B7BA", "#8BB8B7",
                             "#91BAB5", "#95BBB1", "#99BDAD", "#9EBFA9", "#A2C0A5",
                             "#A6C2A0", "#ABC399", "#B0C493", "#B5C68C", "#BAC785",
                             "#BFC87C", "#C5C872", "#CCC968", "#D2CA5D", "#D8CA53",
                             "#DCC847", "#DEC43B", "#DFC02E", "#E1BC21", "#E2B815",
                             "#E3B30E", "#E4AC0C", "#E5A60A", "#E5A007", "#E69905",
                             "#E79305", "#E88D05", "#E98705", "#EA8105", "#EB7B05",
                             "#EC7404", "#ED6D04", "#ED6603", "#EE5E03", "#EE5703",
                             "#EF4B02", "#EF3F01", "#F03301", "#F02700", "#F11B00"])

    row_labels = ["", "A", "B", "C", "D", "E", "F", "G", "H", ""]
    col_labels = [i for i in range(1, 13)]

    # generate coordinates of wells
    coords = [(x, y) for x in range(1, 13) for y in range(1, 9)]

        

    # initialize plot
    fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 8), layout="constrained")
    # set axis limits and ticks
    plt.setp((ax1, ax2), xticks=col_labels, yticklabels=row_labels, xlim=[0, 13], ylim=[0, 9])

    # plot plate 1
    sc = ax1.scatter(*zip(*coords), c=r1_reads[0:96], cmap=zissou, s=250, edgecolors="k", vmin=min(r1_reads), vmax=max(r1_reads))
    ax1.set_title("Plate 1", fontsize=16)
    ax1.invert_yaxis()  # rows labelled top to bottom

    # plot plate 2
    ax2.scatter(*zip(*coords), c=r1_reads[96:192], cmap=zissou, s=250, edgecolors="k")
    ax2.set_title("Plate 2", fontsize=16)
    ax2.invert_yaxis()  # rows labelled top to bottom

    # make color bar
    cbar = fig.colorbar(sc, ax=(ax1, ax2), shrink=0.5)
    cbar.ax.set_title("Reads", fontsize=14)

    plt.savefig(f"{prefix}.png")


def main():
    args = parse_arguments()
    r1_reads_file = getattr(args, "r1_reads_file")
    r1_barcode_file = getattr(args, "r1_barcode_file")
    r1_barcode_set_file = getattr(args, "r1_barcode_set_file")
    prefix = getattr(args, "prefix")

    # read in list of R1 barcodes
    r1_barcodes = [line[0] for line in get_split_lines(r1_barcode_file, delimiter="\t")]

    # map R1 barcodes to their plate coordinates
    # these are the coordinates to be plotted
    r1_coords_dict = get_plate_coords(r1_barcodes)

    # map each observed R1 barcode to its number of reads
    r1_reads_dict = {line[0]: int(line[1]) for line in get_split_lines(r1_reads_file, delimiter="\t", skip=1)}

    # get list of observed reads, in well order
    # these are the values for the color scale
    r1_reads = [r1_reads_dict.get(bc, 0) for bc in r1_barcodes]

    # get nested list of barcodes in each subset
    r1_barcode_subsets = [[line[1:]] for line in get_split_lines(r1_barcode_set_file)]

    # get coordinates of subset splits on plates
    # these are the intercepts of the dashed lines showing plate splits
    plate_divisions = get_plate_divisions(r1_barcode_subsets, r1_coords_dict)

    # make plot
    plate_plot(r1_reads, r1_barcode_subsets, prefix)


if __name__ == "__main__":
    main()
