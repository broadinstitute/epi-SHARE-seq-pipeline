#!/usr/bin/env python
"""
Write output HTML file from list of images and text 

@author Neva Durand (c) 2021
"""
import argparse
import base64
import io
import os.path

def main(output_file_name, image_file_list, log_file_list, input_file_name=None):
    """
    Write to the input file
    Image file list is list of png images
    Log file list is list of text log files

    Separates images by br tag and encodes directly in utf-8
    Log files separated by their title and encoded via pre tag
    """
    # Open output file, write input if exists
    output_file = io.open(output_file_name, 'w', encoding='utf8')
    output_file.write('<!DOCTYPE html><html lang="en"><head><title>Results summary</title><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0"><body>')
    if input_file_name is not None:
        with open(input_file_name) as input_file:
           output_file.write(input_file.read())
    
    with open(image_file_list) as fname:
        images = fname.read().splitlines() 

    # loop through images in image list and encode
    output_file.write('<br>')
    for image in images:
        data = open(image, 'rb').read() # read bytes from file
        data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
        data_base64 = data_base64.decode('utf-8')    # convert bytes to string
        output_file.write('<img width="400" src="data:image/png;base64,' + data_base64 + '" alt=' + os.path.basename(image)+ '><br>') # embed in html

    with open(log_file_list) as fname:
        logs = fname.read().splitlines()

    # loop through log files in log list and write
    for log in logs:
        output_file.write('<h3>')
        output_file.write(os.path.basename(log))
        output_file.write('</h3>')
        output_file.write('<pre>')
        with io.open(log, 'r', encoding='utf8') as log_file:
           text = log_file.read()
        output_file.write(text)
        output_file.write('</pre>')
    output_file.write('</body></html>')
    output_file.close()

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
    group.add_argument('--input_file_name',
                       help='optional file with html text to add at top of file', nargs='?') 
    args = parser.parse_args()
    main(args.output_file_name, args.image_file_list, args.log_file_list, args.input_file_name)

