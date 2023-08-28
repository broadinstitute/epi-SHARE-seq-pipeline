#!/usr/bin/env python
"""
Write output HTML file from list of images and text 

@author Neva Durand (c) 2021
"""
import argparse
import base64
import io
import os.path
import csv

def main(output_file_name, image_file_list, log_file_list, output_csv_name, summary_stats_txt, joint_qc_vals_txt, archr_vals_txt, input_file_name=None):
    """
    Write to the input file
    Image file list is list of png images
    Log file list is list of text log files to link to

    Separates images by br tag and encodes directly in utf-8
    Log files separated by their title and encoded via pre tag
    """
    # Open output file, write input if exists
    output_file = io.open(output_file_name, 'w', encoding='utf8')
    output_file.write('<!DOCTYPE html><html lang="en"><head><title>Results summary</title><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0"><body>')
    if input_file_name is not None:
        with open(input_file_name) as input_file:
           output_file.write(input_file.read())
    
    #csv_output_file = io.open(output_csv_name, 'w', encoding='utf8')
    csv_writer = csv.writer(open(output_csv_name, 'w'))
    
    with open(image_file_list) as fname:
        images = fname.read().splitlines() 

    # loop through images in image list and encode
    output_file.write('<br>')
    print("start of printing images loop")
    for image in images:
        data = open(image, 'rb').read() # read bytes from file
        data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
        data_base64 = data_base64.decode('utf-8')    # convert bytes to string
        name = os.path.basename(image)
        idx = name.index('.') + 1
        name = name[idx:]
        output_file.write('<img id ="' + name + '" width="1000" src="data:image/png;base64,' + data_base64 + '" alt=' + os.path.basename(image)+ '><br>') # embed in html
        #csv_output_file.write(name + ', data:image/png;base64,"' + data_base64 + '\n')
        image_field = '<img id ="' + name + '" width="1000" src="data:image/png;base64,' + data_base64 + '" alt=' + os.path.basename(image)+ '><br>'
        csv_writer.writerow([name, image_field])

    with open(log_file_list) as fname:
        logs = fname.read().splitlines()

    # loop through log files in log list and write
    print("start of printing logs loop")
    for log in logs:
        log_parts = log.split('/')
        log_name = log_parts[-1]
        idx = log_name.index('.') + 1
        log_name = log_name[idx:]
        output_file.write("<div id =" + log_name + ">" + log + "</div>")
        output_file.write("<br>")
        #csv_output_file.write(log_name + ", " + log + '\n')
        csv_writer.writerow([log_name, log])
    output_file.write('</body></html>')
    fname.close()
    
    #write overall summary stats to output csv file
    print("start of printing summary stats loop")
    with open (summary_stats_txt) as fname:
        stats = fname.read().splitlines()
    for stat in stats:
        print("stat is " + stat)
        #csv_output_file.write(stat + '\n')
        name_and_field = stat.split(',')
        name = name_and_field[0]
        stat = name_and_field[1]
        csv_writer.writerow([name, stat])
    fname.close()

    #write stats from qc to output csv file
    print("start of joint qc vals text loop")
    with open (joint_qc_vals_txt) as fname:
        stats = fname.read().splitlines()
    for stat in stats:
        print("stat is " + stat)
        #csv_output_file.write(stat + '\n')
        name_and_field = stat.split(',')
        name = name_and_field[0]
        stat = name_and_field[1]
        csv_writer.writerow([name, stat])
    fname.close()
    
    #write stats from archr to output csv file
    print("print archr vals text loop")
    with open (archr_vals_txt) as fname:
        stats = fname.read().splitlines()
    for stat in stats:
        print("stat is " + stat)
        #csv_output_file.write(stat + '\n')
        name_and_field = stat.split(',')
        name = name_and_field[0]
        stat = name_and_field[1]
        csv_writer.writerow([name, stat])
    fname.close()

    output_file.close()
    #csv_output_file.close()
    csv_writer.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('output_file_name',
                       help='html file to write to')
    group.add_argument('image_file_list', 
                       help='file containing list of image files to paste in HTML file')
    group.add_argument('log_file_list',
                       help='file containing list of text log files to append to end of HTML file')
    group.add_argument('output_csv_name',
                       help='csv file to write values to')
    group.add_argument('summary_stats_txt',
                       help='text file containing names and values of summary stats')
    group.add_argument('joint_qc_vals_txt',
                       help='file containing names and values of stats from joint qc ')
    group.add_argument('archr_vals_txt',
                       help='file containing names and values of stats from archr')
    group.add_argument('--input_file_name',
                       help='optional file with html text to add at top of file', nargs='?') 
    args = parser.parse_args()
    main(args.output_file_name, args.image_file_list, args.log_file_list, args.output_csv_name, args.summary_stats_txt, args.joint_qc_vals_txt, args.archr_vals_txt, args.input_file_name)

