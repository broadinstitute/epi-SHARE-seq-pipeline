#!/usr/bin/env python3

"""
This script takes in a file containing fragment lengths, one per line,
and generates an insert size histogram as a png.
"""

import argparse
import pandas as pd
from collections import defaultdict
from plotnine import *

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot insert size histogram")
    parser.add_argument("fragment_length_file", help="Fragment length file name")
    parser.add_argument("prefix", help="Prefix for plot title")
    parser.add_argument("out_file", help="Name of output png file")
    
    return parser.parse_args()

def get_hist_vals(fragment_length_file):
    """Get dataframe of insert sizes and corresponding counts"""
    with open(fragment_length_file, "r") as f:
        counts = defaultdict(int)
        for line in f:
            fragment_length = int(line.rstrip())
            counts[fragment_length] += 1
            
    df = pd.DataFrame(counts.items(), columns=["insert_size","count"])
    
    return(df)

def label_func(breaks):
    return ["{:.0e}".format(x) for x in breaks]

def plot_hist(df, prefix, out_file):
    plot = (ggplot(df, aes(x="insert_size", y="count")) +
            geom_line(color="red") +
            geom_area(fill="red") +
            labs(title = f"Insert Size Histogram ({prefix})",
                 x = "Insert size",
                 y = "Count") + 
            scale_y_continuous(labels = label_func) +
            theme_classic())
    
    plot.save(filename = out_file, dpi=1000)

def main():
    args = parse_arguments() 
    fragment_length_file = getattr(args, "fragment_length_file")
    prefix = getattr(args, "prefix")
    out_file = getattr(args, "out_file")
    
    df = get_hist_vals(fragment_length_file)
    
    plot_hist(df, prefix, out_file)

if __name__ == "__main__":
    main()    
