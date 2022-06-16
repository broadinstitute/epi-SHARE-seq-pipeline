#!/usr/bin/env python3

"""
Generates a matrix containing read frequency of all barcode combinations
"""

import argparse
import pandas as pd

# Get arguments
desc = "Generates a matrix containing read frequency of all barcode combinations"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-b", help="name of barcode combination file (set of sorted barcodes)")
parser.add_argument("-r", help="name of read frequency file (frequency, barcode)")
parser.add_argument("-o", help="name of output file (R1, R2, R3, PKR, fragments)")

args = parser.parse_args()
observed_barcodes_file = args.b
read_freq_file = args.r
output_file = args.o

# Read in barcodes file
observed_barcodes = pd.read_csv(observed_barcodes_file, header=None, names=["R1","R2","R3","PKR"])
# Get unique values of R1, R2, R3, PKR
R1_set = list(set(observed_barcodes["R1"]))
R2_set = list(set(observed_barcodes["R2"]))
R3_set = list(set(observed_barcodes["R3"]))
PKR_set = list(set(observed_barcodes["PKR"]))

print(f"{len(R1_set)} unique R1 barcodes observed")
print(f"{len(R2_set)} unique R2 barcodes observed")
print(f"{len(R3_set)} unique R3 barcodes observed")
print(f"{len(PKR_set)} unique PKRs observed")

# Generate all combinations of R1, R2, R3, PKR
combos = [(R1,R2,R3,PKR) for R1 in R1_set for R2 in R2_set for R3 in R3_set for PKR in PKR_set]
possible_barcodes = [','.join(x) for x in combos] 

# Make dataframe of R1, R2, R3, PKR, barcodes
bar_freq_df = pd.DataFrame(combos, columns=("R1","R2","R3","PKR"))
bar_freq_df["barcode"] = possible_barcodes

# Read in read frequency file
read_freq = pd.read_csv(read_freq_file, header=None, names=["fragments","barcode"], sep="\t")

# Match possible barcode combinations to observed combinations and get read frequency of matches
bar_freq_df = pd.merge(bar_freq_df, read_freq, how="left", on="barcode")

# Drop barcode column, replace NAs with 0s
bar_freq_df = bar_freq_df.drop("barcode", axis=1)
bar_freq_df = bar_freq_df.fillna(0)

print(bar_freq_df.head())

# Write output file
bar_freq_df.to_csv(output_file, index=False)

