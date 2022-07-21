#!/usr/bin/env python3

"""
Generates a csv containing read frequencies of all possible barcode combinations
"""

import argparse

parser = argparse.ArgumentParser(description="Generates a csv containing read frequencies of all possible barcode combinations")
parser.add_argument("-b", "--barcode_file", help="File containing unique observed barcodes")
parser.add_argument("-r", "--read_freq_file", help="File containing read frequencies of observed barcodes")
parser.add_argument("-o", "--output_file", help="Name of output file")

args = parser.parse_args()

if getattr(args, "barcode_file") is None:
    print("ERROR: Barcode file not provided\n")
    parser.parse_args(["-h"])
    
if getattr(args, "read_freq_file") is None:
    print("ERROR: Read frequency file not provided\n")
    parser.parse_args(["-h"])
    
barcode_file = getattr(args, "barcode_file")
read_freq_file = getattr(args, "read_freq_file")
output_file = getattr(args, "output_file")
    
# use generator to split lines rather than reading into memory
def get_split_lines(file_name, delimiter):
    with open(file_name, "r") as f:
        for line in f:
            yield line.rstrip().split(sep=delimiter)

# create list of all possible barcode combinations
def get_barcode_combos(barcode_file):
    unique_R1 = list({line[0] for line in get_split_lines(barcode_file, delimiter=",")})
    unique_R2 = list({line[1] for line in get_split_lines(barcode_file, delimiter=",")})
    unique_R3 = list({line[2] for line in get_split_lines(barcode_file, delimiter=",")})
    # get unique PKR values only if PKR field exists in barcode file
    unique_PKR = list({line[3] for line in get_split_lines(barcode_file, ",") if len(line) > 3})
    
    if unique_PKR:
        possible_barcodes = [",".join([R1,R2,R3,PKR]) for R1 in unique_R1 for R2 in unique_R2 for R3 in unique_R3 for PKR in unique_PKR]
    else:
        possible_barcodes = [",".join([R1,R2,R3]) for R1 in unique_R1 for R2 in unique_R2 for R3 in unique_R3]
        
    return possible_barcodes    

# create dictionary of read frequencies of all possible barcode combinations               
def get_read_freqs(barcode_file, read_freq_file):
    # make dictionary of observed barcodes and associated read counts
    read_freq_dict = {line[1]:int(line[0]) for line in get_split_lines(read_freq_file, delimiter="\t")}
    
    # initialize dictionary with read counts of 0 for each possible barcode combination
    init_dict = {barcode:0 for barcode in get_barcode_combos(barcode_file)}
    
    # update with counts of observed barcodes
    init_dict.update(read_freq_dict)
    
    return init_dict

read_freqs = get_read_freqs(barcode_file, read_freq_file)

print(f"Writing table of barcode read counts to {output_file}\n")

# write dictionary to csv file    
with open((output_file), "w") as f:
    # write header based on presence/absence of PKR field
    if len(list(read_freqs.keys())[0].split(",")) == 4:        
        f.write("R1,R2,R3,PKR,fragments\n")
    else:
        f.write("R1,R2,R3,fragments\n")
        
    for key in read_freqs.keys():
        f.write("%s,%d\n" % (key, read_freqs[key]))

print("Finished writing csv file\n")