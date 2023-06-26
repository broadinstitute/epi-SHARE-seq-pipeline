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


def main(output_file_name, image_file_list, stats_info, log_file_list, input_file_name=None):
    """
    Write to the input file
    Image file list is list of png images
    Log file list is list of text log files to link to
    Separates images by br tag and encodes directly in utf-8
    Log files separated by their title and encoded via pre tag
    """
    # Open output file, write input if exists
    output_file = io.open(output_file_name, 'w', encoding='utf8')
    
    #function takes a list of pngs, encodes them in base64, and writes them to
    # the output file (all Neva's code)
    def write_pngs(images):
        for image in images:
            data = open(image, 'rb').read() # read bytes from file
            data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
            data_base64 = data_base64.decode('utf-8')    # convert bytes to string
            output_file.write('<img width="1000" src="data:image/png;base64,' + data_base64 + '" alt=' + os.path.basename(image)+ '><br>') # embed in html
    
    #function like write pngs, but ideally works for writing into columns
    def write_pngs_column(images):
        for image in images:
            data = open(image, 'rb').read() # read bytes from file
            data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
            data_base64 = data_base64.decode('utf-8')    # convert bytes to string
            output_file.write('<img style="width:100%;" src="data:image/png;base64,' + data_base64 + '" alt=' + os.path.basename(image)+ '><br>') # embed in html

    
    def write_summary_table(txt, nums, outfile):
        outfile.write("<table>")
        outfile.write("<th>Summary Statistics</th>")
        for index in range(len(txt)):
            if index < len(nums):
                outfile.write("<tr> <td>" + txt[index] + "</td> <td>" + "{:,}".format(float(nums[index])) + " </td> </tr>")
            else:
                outfile.write("<tr> <td>" + txt[index] + "</td> <td>No matching number</td> </tr>")
        outfile.write("</table>")
        
    
    #Sets up the style for the tabs. Also links each tab to the content that
    #should be displayed when the tab is checked
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
    
    output_file.write('<style> table {font-size: 30;} </style>')


    #slighlt adapted from code from following link 
    #https://www.w3schools.com/howto/tryit.asp?filename=tryhow_css_three_columns_unequal
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
    
    
    #this html from the original version, sets up the document
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
    

    # loop through images in image list and encode
    output_file.write('<br>')
    with open(image_file_list) as fname:
        images = fname.read().splitlines() 

    atac_plot_indices = []
    rna_plots_indices = []
    top_level_indices = [0, 1, 17, 24, 26, 29]
    top_level_left_indices = []
    top_level_middle_indicies = []
    top_level_right_indicies = []
    other_indices = []
    
    atac_plots = []
    rna_plots = []
    top_level = [images[i] for i in top_level_indices]
    top_level_left = []
    top_level_middle = []
    top_level_right = []
    other = []
    
    #write the images to the proper tab. The spacing for the output file 
    #write calls is a bit weird because the write_pngs function needs to be 
    #called within a chunk of formatted html
    output_file.write("""<div class="row">
                <div class="column left">""")
    
    # writes summary statistics (eventually should be moved to own tab and 
    # and not hard coded)
    output_file.write("""<center> <span style="font-size: 50;"> 41890 cells </span> <br> <br>
                    <span style="font-size: 30;"> 29656 both 670 RNA 11627 ATAC </span> <br> <br>
                    <span style="font-size: 25;"> RNA:  112,366,191 / 346,626,836 aligned (43% dup) </span> </center>""")
    
    #write images to far left of summary tab 
    write_pngs_column(top_level)
    output_file.write("""</div>
                <div class="column left"> """)
    
    #write images to middle of summary tab
    write_pngs_column(images)
    output_file.write("""</div>
                        <div class="column left">""")
    
    #write images to right of summary tab
    write_pngs_column(images)                    
    output_file.write("""</div>
                    </div>""")
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")
    
    
    #write images to rna tab 
    output_file.write('<div class="tab content2">')
    write_pngs(images)
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")

    #write images to atac plot section
    output_file.write('<div class="tab content3">')
    write_pngs(images)
    output_file.write('<a href="#top">Go to top of page</a>')
    output_file.write("</div>")

    #write to summary stats tab
    stats_names_list = ["Total reads", "Aligned uniquely", "Unaligned", "Unique Reads", "Duplicate Reads", "Percent Duplicates", "NRF=Distinct/Total", "PBC1=OnePair/Distinct", "PBC2=OnePair/TwoPair", "Total reads", "Aligned uniquely", "Aligned multimap", "Unaligned", "Filtered (feature) Reads", "Duplicate Reads", "Percent Duplicates"]
    with open(stats_info) as stats_f:
        stats_list = stats_f.read().split('/n')
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
    group.add_argument('--input_file_name',
                       help='optional file with html text to add at top of file', nargs='?') 

    args = parser.parse_args()
    main(args.output_file_name, args.image_file_list, args.stats_info, args.log_file_list, args.input_file_name)






"""""
bash script to get a file with just the numbers

echo ~{atac_total_reads} "</td></tr>" > sum_stats.txt
echo ~{atac_aligned_uniquely} "</td></tr>" >> sum_stats.txt
echo ~{atac_unaligned} "</td></tr>" >> sum_stats.txt
echo ~{atac_feature_reads} "</td></tr>" >> sum_stats.txt
echo ~{atac_duplicate_reads} "</td></tr>" >> sum_stats.txt
echo ~{atac_percent_duplicates} "</td></tr>" >> sum_stats.txt
echo ~{atac_nrf} "</td></tr>" >> sum_stats.txt
echo ~{atac_pbc1} "</td></tr>" >> sum_stats.txt
echo ~{atac_pbc2} "</td></tr>" >> sum_stats.txt
echo ~{rna_total_reads} "</td></tr>" >> sum_stats.txt
echo ~{rna_aligned_uniquely} "</td></tr>" >> sum_stats.txt
echo ~{rna_aligned_multimap} "</td></tr>" >> sum_stats.txt
echo ~{rna_unaligned} "</td></tr>" >> sum_stats.txt
echo ~{rna_feature_reads} "</td></tr>" >> sum_stats.txt
echo ~{rna_duplicate_reads} "</td></tr>" >> sum_stats.txt
percent=$(( ~{default=0 rna_duplicate_reads}*100/~{default=1 rna_feature_reads} ))
echo $percent "</td></tr></table>" >> sum_stats.txt

names to go with the numbers above
names = ['Total reads', 'Aligned uniquely', 'Unaligned', 'Unique Reads', 'Duplicate Reads', 'Percent Duplicates', 'Distinct/Total', 'OnePair/Distinct', 'OnePair/TwoPair', 'rna switch', 'Total reads', 'Aligned uniquely', 'Aligned multimap', 'Unaligned', 'Filtered', 'Duplicate Reads', 'Percent Duplicates']



#function that writes a table from names and numbers

def write_summary_table(txt, nums, outfile):
    #write the header for the atac section
    outfile.write("<table> <th> ATAC </th>")
    for int in range(len(txt)): 
        #write the header for the rna section after all the atac stuff has 
        #been added
        if int == 6:
            outfile.write("<th> RNA </th>")
        outfile.write("<tr> <td>")
        outfile.write(txt[int])
        outfile.write("</td> <td>")
        #format the numbers to print with commas
        outfile.write(str("{:,}".format(nums[int])))
        outfile.write("</td> </tr>")
        outfile.write("</table>")


#read the text files of numbers and labels into lists, then use those lists
#to make a table
#note: only works if the text files are one label/number per line
label_names_file = open(labels_file, 'r')
#use of splitlines eliminates the trailing endline character
label_names = label_names_file.read().splitlines()
summary_stats_file = open(numbers_file)
summary_stats = summary_stats_file.read().splitlines()
#cast the digits from the text file to ints
summary_stats = [int(num) for num in summary_stats]

write_summary_table(label_names, summary_stats, output_file)

"""

##new thing to push commit 

