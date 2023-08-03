#edited version of the original html print meant to work with the wdl 

#!/usr/bin/env python
"""
Write output HTML file from list of images and text 
@author Neva Durand (c) 2021, edited by Michael Shriver 2023
"""
import argparse
import base64
import io
import os.path
import csv 


def main(output_file_name, image_file_list, stats_info, log_file_list, qc_stats_file, input_file_name=None):
    """
    Write to the input file
    Image file list is list of png images
    Log file list is list of text log files to link to
    Separates images by br tag and encodes directly in utf-8
    Log files separated by their title and encoded via pre tag
    """
    # Open output file, write input if exists
    output_file = io.open(output_file_name, 'w', encoding='utf8')
    
    # function takes a list of pngs, encodes them in base64, and writes them in 
    # image format to an output file at a standard size
    def write_pngs(images):
        for image in images:
            data = open(image, 'rb').read() # read bytes from file
            data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
            data_base64 = data_base64.decode('utf-8')    # convert bytes to string
            output_file.write('<img width="1000" src="data:image/png;base64,' + data_base64 + '" alt=' + os.path.basename(image)+ '><br>') # embed in html
    
    # function takes a list of pngs, encodes them in base64, and writes them in 
    # image format to an output file at the max width possible
    def write_pngs_column(images):
        for image in images:
            data = open(image, 'rb').read() # read bytes from file
            data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
            data_base64 = data_base64.decode('utf-8')    # convert bytes to string
            output_file.write('<img style="width:100%;" src="data:image/png;base64,' + data_base64 + '" alt=' + os.path.basename(image)+ '><br>') # embed in html

    
    #take a string representing a number, and return a string represenation of
    #of that number rounded
    def format_number(txt_num):
        num = float(txt_num)
        if num > 1000000000:
            num = int(num / 1000000)
            num = str(num) + " B"
        elif num > 1000000:
            num = round(num, -6)
            num = int(num / 1000000)
            num = str(num) + " M"
        elif num > 1000:
            num = int(num / 1000)
            num = str(num) + " K"
        elif num < 1: 
            num = str(round(num))
        else: 
            num = str(num)
        return num 
    
    
    # write a table in html to the specified output table with the text data on
    # the right side and the numberic data on the left
    def write_summary_table(txt, nums, outfile):
        outfile.write("<table>")
        outfile.write("<th>Summary Statistics</th>")
        outfile.write("<tr><td colspan=2>ATAC</td></tr>")
        for index in range(len(txt)):
            if index < len(nums):
                #cast string reperesnting number from the text file and print 
                #formated with commas
                outfile.write("<tr> <td>" + txt[index] + "</td> <td>" + format_number(float(nums[index]))+ " </td> </tr>")
                #write rna tab into the table after all of the atac data
                if index == 8:
                    outfile.write("<tr><td colspan=2>RNA</td></tr>")
            else:
                outfile.write("<tr> <td>" + txt[index] + "</td> <td>No matching number</td> </tr>")
        outfile.write("</table>")
        
        
    def write_stats_from_csv(qc_stats_file, outfile, aligned, duplicate):
        file_stats = []
        with open(qc_stats_file, mode='r') as stats_file: 
            file_stats = stats_file.readlines()
        for i in range(len(file_stats)):
            file_stats[i] = ''.join(filter(str.isdigit, file_stats[i]))
  
        total_cells = int(file_stats[0]) + int(file_stats[1]) + int(file_stats[2])
        outfile.write("<center> <span style='font-size: 50;'>" + str(total_cells) + " cells </span> <br> <br>")
        outfile.write("<span style='font-size: 30;'> " + file_stats[0] + " both " + file_stats[1] + " RNA " + file_stats[2] + " ATAC </span> <br> <br>")
        outfile.write("<span style='font-size: 25;'> RNA:   / " +  aligned + " ( " + duplicate + "% dup) </span> </center>")
    
    
    
    # Sets up the style for the tabs. Also links each tab to the content that
    # should be displayed when the tab is checked
    output_file.write("""<style> 
        input {display: none;} 
        input + label {display: inline-block } 
        input ~ .tab {display: none } 
        #tab1:checked ~ .tab.content1,
        #tab2:checked ~ .tab.content2, 
        #tab3:checked ~ .tab.content3,  
        #tab4:checked ~ .tab.content4,
        #tab5:checked ~ .tab.content5{ display: block; } 
        input + label { 
        border: 1px solid #999; 
        background: #EEE; 
        padding: 4px 12px; 
        border-radius: 4px 4px 0 0 ; 
        position: relative; 
        top: 1px; 
        } 
        input:checked + label {
        background: #FFF; 
        border-bottom: 1px solid transparent; 
        } 
        input ~.tab { 
        border-top: 1px solid #999; 
        padding: 12px; 
        } 
        </style>""")
    
    output_file.write('<style> table {font-size: 30; border-spacing; 50px 0} </style>')


    # slightly adapted from code from following link 
    # https://www.w3schools.com/howto/tryit.asp?filename=tryhow_css_three_columns_unequal
    output_file.write(""" <style>
                        * {
                        box-sizing: border-box;
                        }

                        /* Create three unequal columns that floats next to each other */
                        .column {
                        float: left;
                        padding: 10px;
                        }

                        .left {
                        width: 33%;
                        }

                        .right-middle, .right-right {
                        width: 25%;
                        }

                        /* Clear floats after the columns */
                        .row:after {
                        content: "";
                        display: table;
                        clear: both;
                        }
                        
                        </style>
                        """)
    
    
    # write boilerplate to the html file, write information from input file in
    # if input file is provided
    output_file.write('<!DOCTYPE html><html lang="en"><head><title>Results summary</title><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0"><body>')
    if input_file_name is not None:
        with open(input_file_name) as input_file:
           output_file.write(input_file.read())
    
    
    #set up the style and the names for the tabs
    output_file.write("""<input type="radio" name="tabs" id="tab1" checked /> 
        <label for="tab1">Top level Plots</label> 
        <input type="radio" name="tabs" id="tab2" /> 
        <label for="tab2">RNA plots</label> 
        <input type="radio" name="tabs" id="tab3" /"> 
        <label for="tab3">ATAC plots</label>
        <input type="radio" name="tabs" id="tab4" />
        <label for="tab4">Statistics</label>
        <input type="radio" name="tabs" id="tab5" checked />
        <label for="tab5">Log files </label>""")
    
    #everything from here until the first break will be under the first tab
    output_file.write('<div class="tab content1">')
    

    stats_names_list = ["Total reads", "Aligned uniquely", "Unaligned", "Unique Reads", "Duplicate Reads", "Percent Duplicates", "NRF=Distinct/Total", "PBC1=OnePair/Distinct", "PBC2=OnePair/TwoPair", "Total reads", "Aligned uniquely", "Aligned multimap", "Unaligned", "Filtered (feature) Reads", "Duplicate Reads", "Percent Duplicates"]
    with open(stats_info) as stats_f:
        stats_list = stats_f.read().splitlines()
    rna_aligned = stats_list[10]
    rna_duplication_percent = stats_list[15]
    
    # loop through images in image list and encode
    output_file.write('<br>')
    with open(image_file_list) as fname:
        images = fname.read().splitlines() 

    # indicies of plots that should go in each of the tabs
    atac_plot_indices = range(19, len(images))
    rna_plots_indices = range(2, 19)
    top_level_left_indices = [0, 1]
    top_level_middle_indices = [17, 24]
    top_level_right_indices = [26, 29]
    
    # create lists of plots to be written in each tab from indicies
    atac_plots = [images[i] for i in atac_plot_indices]
    rna_plots = [images[i] for i in rna_plots_indices]
    top_level_left = [images[i] for i in top_level_left_indices]
    top_level_middle = [images[i] for i in top_level_middle_indices]
    top_level_right = [images[i] for i in top_level_right_indices]
    
    output_file.write("""<div class="row">
                <div class="column left">""")
    
    # writes summary statistics to left side of summary tab 
    write_stats_from_csv(qc_stats_file, output_file, rna_aligned, rna_duplication_percent)
    
    #write images to left of summary tab 
    write_pngs_column(top_level_left)
    output_file.write("""</div>
                <div class="column left"> """)
    
    #write images to middle of summary tab
    write_pngs_column(top_level_middle)
    output_file.write("""</div>
                        <div class="column left">""")
    
    #write images to right of summary tab
    write_pngs_column(top_level_right)                    
    output_file.write("""</div>
                    </div>""")
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")
    
    
    #write images to rna tab 
    output_file.write('<div class="tab content2">')
    write_pngs(rna_plots)
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")

    #write images to atac plot section
    output_file.write('<div class="tab content3">')
    write_pngs(atac_plots)
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")

    #write to summary stats tab
    output_file.write('<div class="tab content4">') 
    write_summary_table(stats_names_list, stats_list, output_file)
    output_file.write("</div>")
    
    # loop through log files in log list and write
    output_file.write('<div class="tab content5">') 
    with open(log_file_list) as fname:
        logs = fname.read().splitlines()
    for log in logs:
        output_file.write(log)
        output_file.write("<br>")
    output_file.write("</div>")
    output_file.write("</div>")
    output_file.write('</body></html>')
    output_file.close()
    



    
    #need to fix this, going to be something similar to what the images were
    #with open(log_file_list) as fname:
    #comment to push changes

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('output_file_name',
                       help='html file to write to')
    group.add_argument('image_file_list', 
                       help='file containing list of image files to paste in HTML file')
    group.add_argument('stats_info',
                       help='file containing calculated stats')
    group.add_argument('log_file_list',
                       help='file containing list of text log files to append to end of HTML file')
    group.add_argument('qc_stats_file',
                       help='file containing information for top level tab')
    group.add_argument('--input_file_name',
                       help='optional file with html text to add at top of file', nargs='?') 

    args = parser.parse_args()
    main(args.output_file_name, args.image_file_list, args.stats_info, args.log_file_list, args.qc_stats_file, args.input_file_name)