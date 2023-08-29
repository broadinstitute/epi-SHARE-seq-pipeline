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
     # function takes a list of pngs, encodes them in base64, and writes them in 
    # image format to an output file at a standard size
    def write_pngs(images):
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
            #csv_writer.writerow([name, image_field])
            csv_writer.writerow([name, "some image encoding"])
    
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
    
    #I think this will be different, because the file that has the numbers and
    #their labels has them bundled together, so it should be more of a basic
    #line split situation
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
    
    output_file = io.open(output_file_name, 'w', encoding='utf8')
    output_file.write('<!DOCTYPE html><html lang="en"><head><title>Results summary</title><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0"><body>')
    # depricated
    # Open output file, write input if exists
    #if input_file_name is not None:
        #with open(input_file_name) as input_file:
           #output_file.write(input_file.read())
           #output_file.truncate(0)
    
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
    
    #set up the style and the names for the tabs
    output_file.write("""<input type="radio" name="tabs" id="tab1" checked /> 
        <label for="tab1">Top level Plots</label> 
        <input type="radio" name="tabs" id="tab2" /> 
        <label for="tab2">RNA plots</label> 
        <input type="radio" name="tabs" id="tab3" /"> 
        <label for="tab3">ATAC plots</label>
        <input type="radio" name="tabs" id="tab4" />
        <label for="tab4">Statistics</label>
        <input type="radio" name="tabs" id="tab5" />
        <label for="tab5">Log files </label>""")
    
    csv_writer = csv.writer(open(output_csv_name, 'w'))
    
    #everything in this block will be in first tab
    output_file.write('<div class="tab content1">')
    output_file.write('some filler content here')
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")
    
    #write images to second (rna) tab 
    output_file.write('<div class="tab content2">')
    with open(image_file_list) as fname:
        images = fname.read().splitlines() 
        # loop through images in image list and encode
        output_file.write('<br>')
        print("start of printing images loop")
        write_pngs(images)
        output_file.write('<a href="#top">Go to top of page</a>')
        output_file.write("</div>")

    #write images to third (atac) tab
    output_file.write('<div class="tab content3">')
    output_file.write("Plots will go here once things are sorted out")
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")

    #write to fourth (summary stats) tab
    output_file.write('<div class="tab content4">') 
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
            output_file.write(name + " " + stat + '\n')
    fname.close()
    output_file.write("</div>")
    
    # loop through log files in log list and write
    output_file.write('<div class="tab content5">') 
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
    output_file.write("</div>")
    output_file.write("</div>")
    output_file.write('</body></html>')
    output_file.close()     
    
    
    #with open(image_file_list) as fname:
        #images = fname.read().splitlines() 

    # loop through images in image list and encode
    #output_file.write('<br>')
    #print("start of printing images loop")
    

    #with open(log_file_list) as fname:
        #logs = fname.read().splitlines()

    # loop through log files in log list and write
    #print("start of printing logs loop")
    #for log in logs:
        #log_parts = log.split('/')
        #log_name = log_parts[-1]
        #idx = log_name.index('.') + 1
        #log_name = log_name[idx:]
        #output_file.write("<div id =" + log_name + ">" + log + "</div>")
        #output_file.write("<br>")
        #csv_output_file.write(log_name + ", " + log + '\n')
        #csv_writer.writerow([log_name, log])
    #output_file.write('</body></html>')
    #fname.close()
    
    #write overall summary stats to output csv file
    #print("start of printing summary stats loop")
    #with open (summary_stats_txt) as fname:
        #stats = fname.read().splitlines()
    #for stat in stats:
        #print("stat is " + stat)
        #csv_output_file.write(stat + '\n')
        #name_and_field = stat.split(',')
        #name = name_and_field[0]
        #stat = name_and_field[1]
        #csv_writer.writerow([name, stat])
    #fname.close()

    #start csv output only
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
    #fname.close()
    
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
    #fname.close()

    output_file.close()
    #csv_output_file.close()

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

