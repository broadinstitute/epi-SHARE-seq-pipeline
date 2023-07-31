#!/usr/bin/env python3

"""
This script takes in a file containing fragment sizes,
and generates an insert size histogram as a png.
"""

import argparse
import pandas as pd
from collections import defaultdict
from plotnine import *

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot insert size histogram")
    parser.add_argument("fragment_size_file", help="File containing fragment sizes")
    parser.add_argument("prefix", help="Prefix for labelling output file")
    parser.add_argument("out_file", help="Filename for output histogram png")
    
    return parser.parse_args()

def get_hist_vals(fragment_size_file):
    """Get dataframe of histogram values"""
    fragment_size_counts = defaultdict(int)
    with open(fragment_size_file, "r") as f:
        for line in f:
            fragment_size = int(line.rstrip())
            fragment_size_counts[fragment_size] += 1
            
    df = pd.DataFrame(fragment_size_counts.items(), columns=["insert_size","count"])
    
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
    fragment_size_file = getattr(args, "fragment_size_file")
    prefix = getattr(args, "prefix")
    out_file = getattr(args, "out_file")
    
    df = get_hist_vals(fragment_size_file)
    
    plot_hist(df, prefix, out_file)

if __name__ == "__main__":
    main()    
