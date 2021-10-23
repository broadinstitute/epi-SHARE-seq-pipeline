#!/usr/bin/env python

# Author: Sai Ma
# The following program will process the fastq
# handle QC vs full & modify fastq header & split fastqs & add index & trim & split project
# example python3 /mnt/users/sai/Script/Split-seq_Sai/fastq.process.py3.py -a Undetermined_S0_R1_001.fastq.gz -b Undetermined_S0_R2_001.fastq.gz --qc -y /mnt/users/sai/Script/Split-seq_Sai/config_test.yaml

# to do list
## need to add N6 compatability

##### IMPORT MODULES #####
# import necessary for python
import os
import pyfastx
import sys
import time
import xopen
import yaml
from optparse import OptionParser

##### DEFINE FUNCTIONS #####
# Reverse complement
complement = str.maketrans('ATCGN', 'TAGCN')

def _reverse_complement(sequence):
    return sequence.upper().translate(complement)[::-1]

def barcodeSet(barcode, reverse_complement=False):
    bases = "ATCGN"
    barcodeSet = set()
    
    barcodeSet.add(barcode)
    
    if reverse_complement:
        barcode_rc = _reverse_complement(barcode)
        barcodeSet.add(barcode_rc)
    
    for i, nuc in enumerate(barcode):
        if nuc in bases:
            for base in bases:
                if nuc != base:
                    barcodeSet.add((barcode[:i] + base + barcode[i + 1:]))
                    # allow 1 mismatch or 1 bp shift
                    barcodeSet.add((barcode[1:] + base))
                    barcodeSet.add((base + barcode[:-1]))
                    
                    if reverse_complement:
                        barcodeSet.add((barcode_rc[:i] + base + barcode_rc[i + 1:]))
                        # allow 1 mismatch or 1 bp shift
                        barcodeSet.add((barcode_rc[1:] + base))
                        barcodeSet.add((base + barcode_rc[:-1]))
                    
    return barcodeSet
        
P5fwdstr = {'TAGATCGC': 'P1.01',
      'CTCTCTAT': 'P1.02',
      'TATCCTCT': 'P1.03',
      'AGAGTAGA': 'P1.04',
      'GTAAGGAG': 'P1.05',
      'ACTGCATA': 'P1.06',
      'AAGGAGTA': 'P1.07',
      'CTAAGCCT': 'P1.08',
      'TGGAAATC': 'P1.09',
      'AACATGAT': 'P1.10',
      'TGATGAAA': 'P1.11',
      'GTCGGACT': 'P1.12',
      'TTTCTAGC': 'P1.13',
      'TAACCAAG': 'P1.14',
      'GTGTATCG': 'P1.15',
      'TCCATCAA': 'P1.16',
      'TTCGTGCA': 'P1.17',
      'AGGTTGCC': 'P1.18',
      'CCTTATGT': 'P1.19',
      'CAGCAACG': 'P1.20',
      'GGTTCAAT': 'P1.21',
      'ACATTCGT': 'P1.22',
      'GATTCCCA': 'P1.23',
      'CGGACTGC': 'P1.24',
      'AGCCGTTC': 'P1.25',
      'ATTGGGTC': 'P1.26',
      'TGCATACT': 'P1.27',
      'GGGCTTGG': 'P1.28',
      'GACGTGGC': 'P1.29',
      'GCAAATTT': 'P1.30',
      'GCAGCCTC': 'P1.31',
      'TCCGAGTT': 'P1.32',
      'GCATTAAG': 'P1.33',
      'ACGATAAC': 'P1.34',
      'CCTGCGGG': 'P1.35',
      'TGATTGTT': 'P1.36',
      'GGCACGGA': 'P1.37',
      'GATCATTC': 'P1.38',
      'ATGGTCAT': 'P1.39',
      'CGTACCAA': 'P1.40',
      'CCAGTTTA': 'P1.41',
      'ACCGGCCC': 'P1.42',
      'CTAGAAGT': 'P1.43',
      'CGCCAGAT': 'P1.44',
      'TCACATGG': 'P1.45',
      'GAACTCGA': 'P1.46',
      'CCACCGTT': 'P1.47',
      'TAAGTTAC': 'P1.48',
      'GAGACGTG': 'P1.49',
      'TTGCCTAA': 'P1.50',
      'TTAACTTG': 'P1.51',
      'CTTTAACA': 'P1.52',
      'CGTAGACC': 'P1.53',
      'TATTTGCG': 'P1.54',
      'ATCCAGGA': 'P1.55',
      'TGTTCCTG': 'P1.56',
      'ACGCGCAG': 'P1.57',
      'TCTGGCGA': 'P1.58',
      'AATCTACA': 'P1.59',
      'TACTGACC': 'P1.60',
      'CGATAGGG': 'P1.61',
      'ACTTAGAA': 'P1.62',
      'AGAGATCT': 'P1.63',
      'GGTGAAGG': 'P1.64',
      'ATCGAATG': 'P1.65',
      'TCAAGAGC': 'P1.66',
      'GCCCACGT': 'P1.67',
      'TGGGCGGT': 'P1.68',
      'CCCTTGGA': 'P1.69',
      'ATTACCGT': 'P1.70',
      'AGTCCGAG': 'P1.71',
      'ACTTGTTG': 'P1.72',
      'GTAATACA': 'P1.73',
      'GGCGTCTA': 'P1.74',
      'GCGCTGCT': 'P1.75',
      'GTGCCATT': 'P1.76',
      'TAGGTATG': 'P1.77',
      'AACACCTA': 'P1.78',
      'CTCCGAAC': 'P1.79',
      'CAACGGCA': 'P1.80',
      'CAATGTAG': 'P1.81',
      'GGCTACCC': 'P1.82',
      'AAAGTCCG': 'P1.83',
      'TTCCGCGG': 'P1.84',
      'AGGCACTT': 'P1.85',
      'CTTCAGTG': 'P1.86',
      'GCCGGTAG': 'P1.87',
      'TTCAATCC': 'P1.88',
      'CCACACAC': 'P1.89',
      'ATATTATC': 'P1.90',
      'CCGAAGCA': 'P1.91',
      'GTATCGGT': 'P1.92'}

#### OPTIONS ####
# define options
opts = OptionParser()
usage = "usage: %prog [options] [inputs] This will trim adapters"
opts = OptionParser(usage=usage)
opts.add_option("-y", help="<Yaml> yaml file that specifies P5 barcode for each project")
opts.add_option("-a", help="<Read1> Accepts xxx_S1_R1_001.fastq.gz")
opts.add_option("-b", help="<Read2> Accepts xxx_S1_R4_001.fastq.gz")
opts.add_option("--c", help="<Index1> optional xxx_S1_R2_001.fastq.gz", default="NA")
opts.add_option("--d", help="<Index2> optional xxx_S1_R3_001.fastq.gz", default="NA")
#opts.add_option("--outdir", help="output dir")
opts.add_option("--qc", action="store_true", help="QC run with first 3M reads")
opts.add_option("--out", help="Path to the output fastq files")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### INPUTS AND OUTPUTS #####
# name input and outputs
p1_in = options.a
p2_in = options.b
i1_in = options.c
i2_in = options.d
yaml_in = options.y
prefix = options.out
qcreads=99999 # number of reads for QC analysis

 
# name outputs and print to working dir
p1_file = p1_in.split('/')[-1]
p2_file = p2_in.split('/')[-1]

#check for file type and open input file
append = p1_in.split('.')[-1]



# files for reads not matching the P5
dis1 = xopen.xopen('discard/discard.P5.R1.fq.gz', mode='wb')
dis2 = xopen.xopen('discard/discard.P5.R2.fq.gz', mode='wb')
dis3 = xopen.xopen('discard/discard.P5.I1.fq.gz', mode='wb')
dis4 = xopen.xopen('discard/discard.P5.I2.fq.gz', mode='wb')




##### SCRIPT #####
# initialize variables
total_reads = 0
good_reads = 0

# check if reads are indexed
p1_rds = xopen.xopen(p1_in, mode='rb')
p1_line = p1_rds.readline()
seqhead1 = p1_line.decode()

if "+" in seqhead1 or "_" in seqhead1:
    indexed = True
else:
    indexed = False
    if i1_in == "NA":
        sys.exit("warning: one pair of the fastqs are not indexed, but index reads are not supplied")
p1_rds.close()
# end check index


# load yaml 
inFile = open(yaml_in, 'r')
config = yaml.load(inFile, Loader=yaml.UnsafeLoader)
projectNames = set()

for key in config.keys():
    if ('Project' in key):
        projectNames.add(key)
# print(projectNames)

# open files to write in
project = dict()
sampletype = dict()
# N6type = dict()
for proj in projectNames:
    metaData = config[proj]
    outName = metaData['Name']
    outType = metaData['Type']
    # remove existing files
    if os.path.exists("out/mux_" + outName + ".R1.fq.gz"):
        os.remove("out/mux_" + outName + ".R1.fq.gz")
    if os.path.exists("out/mux_" + outName + ".R2.fq.gz"):
        os.remove("out/mux_" + outName + ".R2.fq.gz")

    for primer in metaData['Primer']:
        project[primer] = outName
    for primer in metaData['Primer']:
        sampletype[primer] = outType
            


# Generate i5 barcode sets with 1 mismatch and 1 bp shifting
p5set = dict()
for barcode, name in P5fwdstr.items():
    barcodes = barcodeSet(barcode, reverse_complement=True)
    for bc in barcodes:
        p5set[bc] = name

# Open the connections to all the files that needs to be generated
files_r1 = dict()
files_r2 = dict()
files_i1 = dict()
files_i2 = dict()

for proj in project.values():
    f1 = xopen.xopen("out/mux_" + proj + ".R1.fq.gz", mode='wb')
    f2 = xopen.xopen("out/mux_" + proj + ".R2.fq.gz", mode='wb')
    i1 = xopen.xopen("out/mux_" + proj + ".I1.fq.gz", mode='wb')
    i2 = xopen.xopen("out/mux_" + proj + ".I2.fq.gz", mode='wb')
    files_r1[proj] = f1
    files_r2[proj] = f2
    files_i1[proj] = i1
    files_i2[proj] = i2


start = time.process_time()
open_files = list(map(lambda file: pyfastx.Fastq(file, build_index=False), [p1_in, p2_in, i1_in, i2_in]))

for read1, read2, index1, index2 in zip(*open_files):
    # read1 is a tuple (name, seq, qual)
    # DOES NOT PRINT ANYTHING IN THE HEADER AFTER THE SPACE!
    p5 = index2[1]
    if (p5 in p5set):
            if p5set[p5] in project:
                outName = project[p5set[p5]]
                # Get the file names for the different biosamples
                out_read1 = files_r1[outName]
                out_read2 = files_r2[outName]
                out_index1 = files_i1[outName]
                out_index2 = files_i2[outName]
                
                out_read1.write("@{} 1:N:0:0\n{}\n+\n{}\n".format(read1[0], read1[1], read1[2]).encode())
                out_read2.write("@{} 2:N:0:0\n{}\n+\n{}\n".format(read2[0], read2[1], read2[2]).encode())
                out_index1.write("@{} 1:N:0:0\n{}\n+\n{}\n".format(index1[0], index1[1], index1[2]).encode())
                out_index2.write("@{} 2:N:0:0\n{}\n+\n{}\n".format(index2[0], index2[1], index2[2]).encode())
                
                
                good_reads = good_reads + 1
            else:
                dis1.write("@{} 1:N:0:0\n{}\n+\n{}\n".format(read1[0], read1[1], read1[2]).encode())
                dis2.write("@{} 2:N:0:0\n{}\n+\n{}\n".format(read2[0], read2[1], read2[2]).encode())
                dis3.write("@{} 1:N:0:0\n{}\n+\n{}\n".format(index1[0], index1[1], index1[2]).encode())
                dis4.write("@{} 2:N:0:0\n{}\n+\n{}\n".format(index2[0], index2[1], index2[2]).encode())
    total_reads = total_reads+1

# close files to write the file
map(lambda file: file.close(), files_r1.values())
map(lambda file: file.close(), files_r2.values())
map(lambda file: file.close(), files_i1.values())
map(lambda file: file.close(), files_i2.values())
map(lambda file: file.close(), [dis1, dis2, dis3, dis4])

time = (time.process_time() - start)/60

try:
    print(total_reads)
    print(str(good_reads)+" sucessfully demultiplexed")
    perc = good_reads/(total_reads)*100
    print("%.3g" % perc + "% reads demultiplexed")
    print(f"elapse time:{time}")
except ZeroDivisionError: 
    print("Warning too few reads")
