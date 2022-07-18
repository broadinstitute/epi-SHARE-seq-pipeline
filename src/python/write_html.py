#!/usr/bin/env python
"""
Write output HTML file from list of images and text 

@author Neva Durand (c) 2021
"""
import argparse
import base64

def main(image_files, log_files):
    html = '<br>'

    image_list = image_files.split(",")
    for image in image_list:
        data = open(image, 'rb').read() # read bytes from file
        data_base64 = base64.b64encode(data)  # encode to base64 (bytes)
        data_base64 = data_base64.decode()    # convert bytes to string
        html = html + '<img src="data:image/jpeg;base64,' + data_base64 + '"><br>' # embed in html
    f = open('output.html', 'a+')
    f.write(html)

    f.write('<pre>')
    log_list = log_files.split(",")
    for log in log_list:
        f.write('<h3>')
        f.write(log)
        f.write('</h3>')
        f.write('<pre>')
	f1 = open(log, 'r')
        f.write(f1.read())
        f.write('</pre>)
    f.write('</body></html>')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('image_files', 
                       help='comma separated list of image files to paste in HTML file')
    group.add_argument('log_files',
                       help='comma separated list of text log files to append to end of HTML file')
    args = parser.parse_args()
    main(args.image_files, args.log_files)

