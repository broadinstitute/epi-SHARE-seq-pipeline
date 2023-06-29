#!/usr/bin/env python
"""
Write output CVS file with that has a label and data for important numbers,
files, and logs generated by the pipeline
@author Neva Durand (c) 2021, edited by Michael Shriver 2023
"""

import argparse
import io
import pandas as pd 

def main(output_file_name, names_list, numeric_list, image_list, log_list):

    #output_file = io.open(output_file_name, 'w', encoding='utf8')
    
    # generate an array of names for data
    names_array = []
    if names_list is not None:
        with open(names_list, "r") as names_f: 
            names_array = names_f.read().splitlines()
    names_array_clean = []
    for name in names_array:
        names_array_clean.append((name.split('/')[-1]))    
    
    # generate an array of the numeric data (represented by strings)
    if numeric_list is not None:
        with open(numeric_list, "r") as numeric_f: 
            numeric_array = numeric_f.read().split('/n')
    
    # generate an array of image data (represented by strings)
    if image_list is not None:
        with open(image_list, "r") as image_f: 
            image_array = image_f.read().splitlines()

    # generate an array of logs files
    if log_list is not None:
        with open(log_list, "r") as log_f: 
            log_array = log_f.read().splitlines()

    # combine the arrays that represent values into one array that lines up
    # with the array of names of things they represent
    csv_values = numeric_array + image_array + log_array
    
    # make a csv file with labels formated label, data
    csv_data_dict = {'Name': names_array_clean, 'Data': csv_values}
    csv_df = pd.DataFrame(csv_data_dict)
    csv_df.to_csv(output_file_name)







if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_argument_group()
    group.add_argument('output_file_name',
                       help='html file to write to')
    group.add_argument('names_list', 
                       help='file containing list of names of relevant data from a pipeline run')
    group.add_argument('numeric_list',
                       help='file containing list of numbers produced by the pipeline')
    group.add_argument('image_list',
                       help='file containing list of png strings?') 
    
    group.add_argument('log_list',
                       help='file containing the log information for a run')
    args = parser.parse_args()
    main(args.output_file_name, args.names_list, args.numeric_list, args.image_list, args.log_list)