#!/usr/bin/env python3

import csv
import os
import sys
import yaml
from optparse import OptionParser
#### OPTIONS ####
# define options
opts = OptionParser()
usage = "usage: %prog [options] [inputs] This will trim adapters"
opts = OptionParser(usage=usage)
opts.add_option("-y", help="<Yaml> yaml configuration file path")
opts.add_option("-p", help="<Paths> Comma-separated string containing the path to the fastqs")
opts.add_option("-n", help="<Name> Name to give to the table")
options, arguments = opts.parse_args()

##### INPUTS AND OUTPUTS #####
# name input and outputs
yaml_in = options.y
paths_in = options.p
table_name = options.n

paths_arr = paths_in.strip().split(',')

idx = 0
terra_table_name = f"entity:{table_name}_Id"
fieldnames = [terra_table_name, 'Names', 'Type', 'Primer', 'fastq_R1', 'fastq_R2']

with open('results_load_table.tsv', 'w', newline='') as f_output:
    #csv_output = csv.DictWriter(f_output, fieldnames=fieldnames)
    #csv_output.writeheader()
    f_output.write("\t".join(fieldnames))
    f_output.write("\n")

    with open(yaml_in) as f_input:
        data = yaml.safe_load(f_input)
        for entry in data:
            row = {}
            row[terra_table_name] = entry
            print(data[entry])
            for get in ['Name', 'Type', 'Primer']:
                try:
                    if get == "Name":
                        row['Names'] = data[entry][get]
                    elif get == "Primer":
                        row[get] = ",".join(data[entry][get])
                    else:
                        row[get] = data[entry][get]
                except KeyError as e:
                    pass
            row['fastq_R1'] = paths_arr[idx] 
            row['fastq_R2'] = paths_arr[idx+1]
            f_output.write("\t".join(row.values()))
            f_output.write("\n")
            #csv_output.writerow(row)
            idx += 2


