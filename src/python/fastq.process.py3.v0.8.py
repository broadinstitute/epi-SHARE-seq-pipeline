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
import re
import sys
import bz2
import gzip
import string
import Levenshtein
import json
import yaml
from Bio import SeqIO
from Bio import AlignIO
from optparse import OptionParser
import time
import io

##### DEFINE FUNCTIONS #####
# Reverse complement
complement = str.maketrans('ATCGN', 'TAGCN')

def reverse_complement(sequence):
    return sequence.decode("utf-8").upper().translate(complement)[::-1]

# Align with mismatch, find first and move on, assumes only one
def fuzz_align(s_seq,l_seq,mismatch):
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq)
        if dist <= mismatch:  # find first then break
            return i, dist
            break

def barcodeSet(barcode):
    bases = "ATCGN"
    barcodeSet = set()
    barcodeSet.add(barcode)
    for i, c in enumerate(barcode):
        if c in bases:
            for base in bases:
                if c != base:
                    barcodeSet.add((barcode[:i] + base + barcode[i + 1:]))
                    ## allow 1 mismatch or 1 bp shift
                    # barcodeSet.add((barcode[1:] + base))
                    # barcodeSet.add((base + barcode[:-1]))
    return barcodeSet
        
## next-seq, nova-seq 1.5 P5 barcode
P5revcomp = {'GCGATCTA': 'P1.01',
      'ATAGAGAG': 'P1.02',
      'AGAGGATA': 'P1.03',
      'TCTACTCT': 'P1.04',
      'CTCCTTAC': 'P1.05',
      'TATGCAGT': 'P1.06',
      'TACTCCTT': 'P1.07',
      'AGGCTTAG': 'P1.08',
      'GATTTCCA': 'P1.09',
      'ATCATGTT': 'P1.10',
      'TTTCATCA': 'P1.11',
      'AGTCCGAC': 'P1.12',
      'GCTAGAAA': 'P1.13',
      'CTTGGTTA': 'P1.14',
      'CGATACAC': 'P1.15',
      'TTGATGGA': 'P1.16',
      'TGCACGAA': 'P1.17',
      'GGCAACCT': 'P1.18',
      'ACATAAGG': 'P1.19',
      'CGTTGCTG': 'P1.20',
      'ATTGAACC': 'P1.21',
      'ACGAATGT': 'P1.22',
      'TGGGAATC': 'P1.23',
      'GCAGTCCG': 'P1.24',
      'GAACGGCT': 'P1.25',
      'GACCCAAT': 'P1.26',
      'AGTATGCA': 'P1.27',
      'CCAAGCCC': 'P1.28',
      'GCCACGTC': 'P1.29',
      'AAATTTGC': 'P1.30',
      'GAGGCTGC': 'P1.31',
      'AACTCGGA': 'P1.32',
      'CTTAATGC': 'P1.33',
      'GTTATCGT': 'P1.34',
      'CCCGCAGG': 'P1.35',
      'AACAATCA': 'P1.36',
      'TCCGTGCC': 'P1.37',
      'GAATGATC': 'P1.38',
      'ATGACCAT': 'P1.39',
      'TTGGTACG': 'P1.40',
      'TAAACTGG': 'P1.41',
      'GGGCCGGT': 'P1.42',
      'ACTTCTAG': 'P1.43',
      'ATCTGGCG': 'P1.44',
      'CCATGTGA': 'P1.45',
      'TCGAGTTC': 'P1.46',
      'AACGGTGG': 'P1.47',
      'GTAACTTA': 'P1.48',
      'CACGTCTC': 'P1.49',
      'TTAGGCAA': 'P1.50',
      'CAAGTTAA': 'P1.51',
      'TGTTAAAG': 'P1.52',
      'GGTCTACG': 'P1.53',
      'CGCAAATA': 'P1.54',
      'TCCTGGAT': 'P1.55',
      'CAGGAACA': 'P1.56',
      'CTGCGCGT': 'P1.57',
      'TCGCCAGA': 'P1.58',
      'TGTAGATT': 'P1.59',
      'GGTCAGTA': 'P1.60',
      'CCCTATCG': 'P1.61',
      'TTCTAAGT': 'P1.62',
      'AGATCTCT': 'P1.63',
      'CCTTCACC': 'P1.64',
      'CATTCGAT': 'P1.65',
      'GCTCTTGA': 'P1.66',
      'ACGTGGGC': 'P1.67',
      'ACCGCCCA': 'P1.68',
      'TCCAAGGG': 'P1.69',
      'ACGGTAAT': 'P1.70',
      'CTCGGACT': 'P1.71',
      'CAACAAGT': 'P1.72',
      'TGTATTAC': 'P1.73',
      'TAGACGCC': 'P1.74',
      'AGCAGCGC': 'P1.75',
      'AATGGCAC': 'P1.76',
      'CATACCTA': 'P1.77',
      'TAGGTGTT': 'P1.78',
      'GTTCGGAG': 'P1.79',
      'TGCCGTTG': 'P1.80',
      'CTACATTG': 'P1.81',
      'GGGTAGCC': 'P1.82',
      'CGGACTTT': 'P1.83',
      'CCGCGGAA': 'P1.84',
      'AAGTGCCT': 'P1.85',
      'CACTGAAG': 'P1.86',
      'CTACCGGC': 'P1.87',
      'GGATTGAA': 'P1.88',
      'GTGTGTGG': 'P1.89',
      'GATAATAT': 'P1.90',
      'TGCTTCGG': 'P1.91',
      'ACCGATAC': 'P1.92'}

## mi-seq, nova-seq 1.0 P5 barcode
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
R1 = {'ATCACGTT': 'R1.001',
      'CGATGTTT': 'R1.002',
      'TTAGGCAT': 'R1.003',
      'TGACCACT': 'R1.004',
      'ACAGTGGT': 'R1.005',
      'GCCAATGT': 'R1.006',
      'CAGATCTG': 'R1.007',
      'ACTTGATG': 'R1.008',
      'GATCAGCG': 'R1.009',
      'TAGCTTGT': 'R1.010',
      'GGCTACAG': 'R1.011',
      'CTTGTACT': 'R1.012',
      'TGGTTGTT': 'R1.013',
      'TCTCGGTT': 'R1.014',
      'TAAGCGTT': 'R1.015',
      'TCCGTCTT': 'R1.016',
      'TGTACCTT': 'R1.017',
      'TTCTGTGT': 'R1.018',
      'TCTGCTGT': 'R1.019',
      'TTGGAGGT': 'R1.020',
      'TCGAGCGT': 'R1.021',
      'TGATACGT': 'R1.022',
      'TGCATAGT': 'R1.023',
      'TTGACTCT': 'R1.024',
      'TGCGATCT': 'R1.025',
      'TTCCTGCT': 'R1.026',
      'TAGTGACT': 'R1.027',
      'TACAGGAT': 'R1.028',
      'TCCTCAAT': 'R1.029',
      'TGTGGTTG': 'R1.030',
      'TACTAGTC': 'R1.031',
      'TTCCATTG': 'R1.032',
      'TCGAAGTG': 'R1.033',
      'TAACGCTG': 'R1.034',
      'TTGGTATG': 'R1.035',
      'TGAACTGG': 'R1.036',
      'TACTTCGG': 'R1.037',
      'TCTCACGG': 'R1.038',
      'TCAGGAGG': 'R1.039',
      'TAAGTTCG': 'R1.040',
      'TCCAGTCG': 'R1.041',
      'TGTATGCG': 'R1.042',
      'TCATTGAG': 'R1.043',
      'TGGCTCAG': 'R1.044',
      'TATGCCAG': 'R1.045',
      'TCAGATTC': 'R1.046',
      'TAGTCTTG': 'R1.047',
      'TTCAGCTC': 'R1.048',
      'TGTCTATC': 'R1.049',
      'TATGTGGC': 'R1.050',
      'TTACTCGC': 'R1.051',
      'TCGTTAGC': 'R1.052',
      'TACCGAGC': 'R1.053',
      'TGTTCTCC': 'R1.054',
      'TTCGCACC': 'R1.055',
      'TTGCGTAC': 'R1.056',
      'TCTACGAC': 'R1.057',
      'TGACAGAC': 'R1.058',
      'TAGAACAC': 'R1.059',
      'TCATCCTA': 'R1.060',
      'TGCTGATA': 'R1.061',
      'TAGACGGA': 'R1.062',
      'TGTGAAGA': 'R1.063',
      'TCTCTTCA': 'R1.064',
      'TTGTTCCA': 'R1.065',
      'TGAAGCCA': 'R1.066',
      'TACCACCA': 'R1.067',
      'TGCGTGAA': 'R1.068',
      'GGTGAGTT': 'R1.069',
      'GATCTCTT': 'R1.070',
      'GTGTCCTT': 'R1.071',
      'GACGGATT': 'R1.072',
      'GCAACATT': 'R1.073',
      'GGTCGTGT': 'R1.074',
      'GAATCTGT': 'R1.075',
      'GTACATCT': 'R1.076',
      'GAGGTGCT': 'R1.077',
      'GCATGGCT': 'R1.078',
      'GTTAGCCT': 'R1.079',
      'GTCGCTAT': 'R1.080',
      'GGAATGAT': 'R1.081',
      'GAGCCAAT': 'R1.082',
      'GCTCCTTG': 'R1.083',
      'GTAAGGTG': 'R1.084',
      'GAGGATGG': 'R1.085',
      'GTTGTCGG': 'R1.086',
      'GGATTAGG': 'R1.087',
      'GATAGAGG': 'R1.088',
      'GTGTGTCG': 'R1.089',
      'GCAATCCG': 'R1.090',
      'GACCTTAG': 'R1.091',
      'GCCTGTTC': 'R1.092',
      'GCACTGTC': 'R1.093',
      'GCTAACTC': 'R1.094',
      'GATTCATC': 'R1.095',
      'GTCTTGGC': 'R1.096',
      'ACTGGTAG': 'R1.097',
      'GCTCCAAC': 'R1.098',
      'GCGTAAGA': 'R1.099',
      'TGACCATC': 'R1.100',
      'GGATTCTC': 'R1.101',
      'TCGTTGCG': 'R1.102',
      'GCCAGACT': 'R1.103',
      'TAGGTCCG': 'R1.104',
      'CCGCCTCC': 'R1.105',
      'AGTAGATC': 'R1.106',
      'AGTCTTAA': 'R1.107',
      'GAATGAAT': 'R1.108',
      'CTGAAGTC': 'R1.109',
      'TTATAGCT': 'R1.110',
      'CTTCTCTT': 'R1.111',
      'CCTATTGA': 'R1.112',
      'GCGGCGCA': 'R1.113',
      'AGGCAGTA': 'R1.114',
      'GTTAATTC': 'R1.115',
      'GGAGAATG': 'R1.116',
      'TCTGCCGG': 'R1.117',
      'GCCTCGGT': 'R1.118',
      'CTAATATG': 'R1.119',
      'TGGCTGAC': 'R1.120',
      'ATTGCATA': 'R1.121',
      'AGTCGCGT': 'R1.122',
      'GTTCCGTA': 'R1.123',
      'ATCCGGAT': 'R1.124',
      'TATTACTC': 'R1.125',
      'ACGGACCT': 'R1.126',
      'GTAAGGCC': 'R1.127',
      'GAGACGTT': 'R1.128',
      'TTGGAACT': 'R1.129',
      'AAGGTTAA': 'R1.130',
      'TTAGCGGT': 'R1.131',
      'TATCAAGG': 'R1.132',
      'TGGTATCA': 'R1.133',
      'CTCGATCT': 'R1.134',
      'CTAGCTTC': 'R1.135',
      'CGAGGCCG': 'R1.136',
      'AACCTTCT': 'R1.137',
      'TCGGAGTA': 'R1.138',
      'GCAAGTTA': 'R1.139',
      'GAGCTCAG': 'R1.140',
      'CAGGAATC': 'R1.141',
      'GGTTGGCG': 'R1.142',
      'TGCTCCTT': 'R1.143',
      'CCGGCGTT': 'R1.144',
      'CAAGCATA': 'R1.145',
      'TACGCAGG': 'R1.146',
      'GGCGTCTG': 'R1.147',
      'ATGACGGC': 'R1.148',
      'AGACTACT': 'R1.149',
      'AATCTCGC': 'R1.150',
      'GGCGGATA': 'R1.151',
      'GCGCCTAG': 'R1.152',
      'TGCTTAGC': 'R1.153',
      'TTAGTATA': 'R1.154',
      'TCTAGCAA': 'R1.155',
      'TCCAATGG': 'R1.156',
      'CCAACCGC': 'R1.157',
      'GGTAACTA': 'R1.158',
      'GATGATTG': 'R1.159',
      'ATTATCGA': 'R1.160',
      'TAGATAAT': 'R1.161',
      'CTTACGTT': 'R1.162',
      'TACTATTA': 'R1.163',
      'AGAAGCGA': 'R1.164',
      'ATCTAGTT': 'R1.165',
      'AAGTTCCA': 'R1.166',
      'GAACTAAC': 'R1.167',
      'CCGTCGAA': 'R1.168',
      'TGATTGAA': 'R1.169',
      'TTCTAATC': 'R1.170',
      'GATAGCCG': 'R1.171',
      'TTGACTTC': 'R1.172',
      'GTCGAGGT': 'R1.173',
      'GCTATCCT': 'R1.174',
      'ATAGAGAG': 'R1.175',
      'TTCTTGAC': 'R1.176',
      'GTTAGGAG': 'R1.177',
      'AGTCATGG': 'R1.178',
      'AGACCAGG': 'R1.179',
      'TGCCGGTA': 'R1.180',
      'TAACTGGC': 'R1.181',
      'GTCGGCAA': 'R1.182',
      'CAACGGTC': 'R1.183',
      'ATGCTAAC': 'R1.184',
      'TGGAGGCG': 'R1.185',
      'CCTTGGCT': 'R1.186',
      'GGATGGTA': 'R1.187',
      'GAGGAGAG': 'R1.188',
      'CGGATAAG': 'R1.189',
      'TGCAGAGA': 'R1.190',
      'TTGGCCAG': 'R1.191',
      'TAACCAGT': 'R1.192'}
R2 = {'ATCACGTT': 'R2.001',
      'CGATGTTT': 'R2.002',
      'TTAGGCAT': 'R2.003',
      'TGACCACT': 'R2.004',
      'ACAGTGGT': 'R2.005',
      'GCCAATGT': 'R2.006',
      'CAGATCTG': 'R2.007',
      'ACTTGATG': 'R2.008',
      'GATCAGCG': 'R2.009',
      'TAGCTTGT': 'R2.010',
      'GGCTACAG': 'R2.011',
      'CTTGTACT': 'R2.012',
      'TGGTTGTT': 'R2.013',
      'TCTCGGTT': 'R2.014',
      'TAAGCGTT': 'R2.015',
      'TCCGTCTT': 'R2.016',
      'TGTACCTT': 'R2.017',
      'TTCTGTGT': 'R2.018',
      'TCTGCTGT': 'R2.019',
      'TTGGAGGT': 'R2.020',
      'TCGAGCGT': 'R2.021',
      'TGATACGT': 'R2.022',
      'TGCATAGT': 'R2.023',
      'TTGACTCT': 'R2.024',
      'TGCGATCT': 'R2.025',
      'TTCCTGCT': 'R2.026',
      'TAGTGACT': 'R2.027',
      'TACAGGAT': 'R2.028',
      'TCCTCAAT': 'R2.029',
      'TGTGGTTG': 'R2.030',
      'TACTAGTC': 'R2.031',
      'TTCCATTG': 'R2.032',
      'TCGAAGTG': 'R2.033',
      'TAACGCTG': 'R2.034',
      'TTGGTATG': 'R2.035',
      'TGAACTGG': 'R2.036',
      'TACTTCGG': 'R2.037',
      'TCTCACGG': 'R2.038',
      'TCAGGAGG': 'R2.039',
      'TAAGTTCG': 'R2.040',
      'TCCAGTCG': 'R2.041',
      'TGTATGCG': 'R2.042',
      'TCATTGAG': 'R2.043',
      'TGGCTCAG': 'R2.044',
      'TATGCCAG': 'R2.045',
      'TCAGATTC': 'R2.046',
      'TAGTCTTG': 'R2.047',
      'TTCAGCTC': 'R2.048',
      'TGTCTATC': 'R2.049',
      'TATGTGGC': 'R2.050',
      'TTACTCGC': 'R2.051',
      'TCGTTAGC': 'R2.052',
      'TACCGAGC': 'R2.053',
      'TGTTCTCC': 'R2.054',
      'TTCGCACC': 'R2.055',
      'TTGCGTAC': 'R2.056',
      'TCTACGAC': 'R2.057',
      'TGACAGAC': 'R2.058',
      'TAGAACAC': 'R2.059',
      'TCATCCTA': 'R2.060',
      'TGCTGATA': 'R2.061',
      'TAGACGGA': 'R2.062',
      'TGTGAAGA': 'R2.063',
      'TCTCTTCA': 'R2.064',
      'TTGTTCCA': 'R2.065',
      'TGAAGCCA': 'R2.066',
      'TACCACCA': 'R2.067',
      'TGCGTGAA': 'R2.068',
      'GGTGAGTT': 'R2.069',
      'GATCTCTT': 'R2.070',
      'GTGTCCTT': 'R2.071',
      'GACGGATT': 'R2.072',
      'GCAACATT': 'R2.073',
      'GGTCGTGT': 'R2.074',
      'GAATCTGT': 'R2.075',
      'GTACATCT': 'R2.076',
      'GAGGTGCT': 'R2.077',
      'GCATGGCT': 'R2.078',
      'GTTAGCCT': 'R2.079',
      'GTCGCTAT': 'R2.080',
      'GGAATGAT': 'R2.081',
      'GAGCCAAT': 'R2.082',
      'GCTCCTTG': 'R2.083',
      'GTAAGGTG': 'R2.084',
      'GAGGATGG': 'R2.085',
      'GTTGTCGG': 'R2.086',
      'GGATTAGG': 'R2.087',
      'GATAGAGG': 'R2.088',
      'GTGTGTCG': 'R2.089',
      'GCAATCCG': 'R2.090',
      'GACCTTAG': 'R2.091',
      'GCCTGTTC': 'R2.092',
      'GCACTGTC': 'R2.093',
      'GCTAACTC': 'R2.094',
      'GATTCATC': 'R2.095',
      'GTCTTGGC': 'R2.096',
      'ACTGGTAG': 'R2.097',
      'GCTCCAAC': 'R2.098',
      'GCGTAAGA': 'R2.099',
      'TGACCATC': 'R2.100',
      'GGATTCTC': 'R2.101',
      'TCGTTGCG': 'R2.102',
      'GCCAGACT': 'R2.103',
      'TAGGTCCG': 'R2.104',
      'CCGCCTCC': 'R2.105',
      'AGTAGATC': 'R2.106',
      'AGTCTTAA': 'R2.107',
      'GAATGAAT': 'R2.108',
      'CTGAAGTC': 'R2.109',
      'TTATAGCT': 'R2.110',
      'CTTCTCTT': 'R2.111',
      'CCTATTGA': 'R2.112',
      'GCGGCGCA': 'R2.113',
      'AGGCAGTA': 'R2.114',
      'GTTAATTC': 'R2.115',
      'GGAGAATG': 'R2.116',
      'TCTGCCGG': 'R2.117',
      'GCCTCGGT': 'R2.118',
      'CTAATATG': 'R2.119',
      'TGGCTGAC': 'R2.120',
      'ATTGCATA': 'R2.121',
      'AGTCGCGT': 'R2.122',
      'GTTCCGTA': 'R2.123',
      'ATCCGGAT': 'R2.124',
      'TATTACTC': 'R2.125',
      'ACGGACCT': 'R2.126',
      'GTAAGGCC': 'R2.127',
      'GAGACGTT': 'R2.128',
      'TTGGAACT': 'R2.129',
      'AAGGTTAA': 'R2.130',
      'TTAGCGGT': 'R2.131',
      'TATCAAGG': 'R2.132',
      'TGGTATCA': 'R2.133',
      'CTCGATCT': 'R2.134',
      'CTAGCTTC': 'R2.135',
      'CGAGGCCG': 'R2.136',
      'AACCTTCT': 'R2.137',
      'TCGGAGTA': 'R2.138',
      'GCAAGTTA': 'R2.139',
      'GAGCTCAG': 'R2.140',
      'CAGGAATC': 'R2.141',
      'GGTTGGCG': 'R2.142',
      'TGCTCCTT': 'R2.143',
      'CCGGCGTT': 'R2.144',
      'CAAGCATA': 'R2.145',
      'TACGCAGG': 'R2.146',
      'GGCGTCTG': 'R2.147',
      'ATGACGGC': 'R2.148',
      'AGACTACT': 'R2.149',
      'AATCTCGC': 'R2.150',
      'GGCGGATA': 'R2.151',
      'GCGCCTAG': 'R2.152',
      'TGCTTAGC': 'R2.153',
      'TTAGTATA': 'R2.154',
      'TCTAGCAA': 'R2.155',
      'TCCAATGG': 'R2.156',
      'CCAACCGC': 'R2.157',
      'GGTAACTA': 'R2.158',
      'GATGATTG': 'R2.159',
      'ATTATCGA': 'R2.160',
      'TAGATAAT': 'R2.161',
      'CTTACGTT': 'R2.162',
      'TACTATTA': 'R2.163',
      'AGAAGCGA': 'R2.164',
      'ATCTAGTT': 'R2.165',
      'AAGTTCCA': 'R2.166',
      'GAACTAAC': 'R2.167',
      'CCGTCGAA': 'R2.168',
      'TGATTGAA': 'R2.169',
      'TTCTAATC': 'R2.170',
      'GATAGCCG': 'R2.171',
      'TTGACTTC': 'R2.172',
      'GTCGAGGT': 'R2.173',
      'GCTATCCT': 'R2.174',
      'ATAGAGAG': 'R2.175',
      'TTCTTGAC': 'R2.176',
      'GTTAGGAG': 'R2.177',
      'AGTCATGG': 'R2.178',
      'AGACCAGG': 'R2.179',
      'TGCCGGTA': 'R2.180',
      'TAACTGGC': 'R2.181',
      'GTCGGCAA': 'R2.182',
      'CAACGGTC': 'R2.183',
      'ATGCTAAC': 'R2.184',
      'TGGAGGCG': 'R2.185',
      'CCTTGGCT': 'R2.186',
      'GGATGGTA': 'R2.187',
      'GAGGAGAG': 'R2.188',
      'CGGATAAG': 'R2.189',
      'TGCAGAGA': 'R2.190',
      'TTGGCCAG': 'R2.191',
      'TAACCAGT': 'R2.192'}
R3 = {'ATCACGTT': 'R3.001',
      'CGATGTTT': 'R3.002',
      'TTAGGCAT': 'R3.003',
      'TGACCACT': 'R3.004',
      'ACAGTGGT': 'R3.005',
      'GCCAATGT': 'R3.006',
      'CAGATCTG': 'R3.007',
      'ACTTGATG': 'R3.008',
      'GATCAGCG': 'R3.009',
      'TAGCTTGT': 'R3.010',
      'GGCTACAG': 'R3.011',
      'CTTGTACT': 'R3.012',
      'TGGTTGTT': 'R3.013',
      'TCTCGGTT': 'R3.014',
      'TAAGCGTT': 'R3.015',
      'TCCGTCTT': 'R3.016',
      'TGTACCTT': 'R3.017',
      'TTCTGTGT': 'R3.018',
      'TCTGCTGT': 'R3.019',
      'TTGGAGGT': 'R3.020',
      'TCGAGCGT': 'R3.021',
      'TGATACGT': 'R3.022',
      'TGCATAGT': 'R3.023',
      'TTGACTCT': 'R3.024',
      'TGCGATCT': 'R3.025',
      'TTCCTGCT': 'R3.026',
      'TAGTGACT': 'R3.027',
      'TACAGGAT': 'R3.028',
      'TCCTCAAT': 'R3.029',
      'TGTGGTTG': 'R3.030',
      'TACTAGTC': 'R3.031',
      'TTCCATTG': 'R3.032',
      'TCGAAGTG': 'R3.033',
      'TAACGCTG': 'R3.034',
      'TTGGTATG': 'R3.035',
      'TGAACTGG': 'R3.036',
      'TACTTCGG': 'R3.037',
      'TCTCACGG': 'R3.038',
      'TCAGGAGG': 'R3.039',
      'TAAGTTCG': 'R3.040',
      'TCCAGTCG': 'R3.041',
      'TGTATGCG': 'R3.042',
      'TCATTGAG': 'R3.043',
      'TGGCTCAG': 'R3.044',
      'TATGCCAG': 'R3.045',
      'TCAGATTC': 'R3.046',
      'TAGTCTTG': 'R3.047',
      'TTCAGCTC': 'R3.048',
      'TGTCTATC': 'R3.049',
      'TATGTGGC': 'R3.050',
      'TTACTCGC': 'R3.051',
      'TCGTTAGC': 'R3.052',
      'TACCGAGC': 'R3.053',
      'TGTTCTCC': 'R3.054',
      'TTCGCACC': 'R3.055',
      'TTGCGTAC': 'R3.056',
      'TCTACGAC': 'R3.057',
      'TGACAGAC': 'R3.058',
      'TAGAACAC': 'R3.059',
      'TCATCCTA': 'R3.060',
      'TGCTGATA': 'R3.061',
      'TAGACGGA': 'R3.062',
      'TGTGAAGA': 'R3.063',
      'TCTCTTCA': 'R3.064',
      'TTGTTCCA': 'R3.065',
      'TGAAGCCA': 'R3.066',
      'TACCACCA': 'R3.067',
      'TGCGTGAA': 'R3.068',
      'GGTGAGTT': 'R3.069',
      'GATCTCTT': 'R3.070',
      'GTGTCCTT': 'R3.071',
      'GACGGATT': 'R3.072',
      'GCAACATT': 'R3.073',
      'GGTCGTGT': 'R3.074',
      'GAATCTGT': 'R3.075',
      'GTACATCT': 'R3.076',
      'GAGGTGCT': 'R3.077',
      'GCATGGCT': 'R3.078',
      'GTTAGCCT': 'R3.079',
      'GTCGCTAT': 'R3.080',
      'GGAATGAT': 'R3.081',
      'GAGCCAAT': 'R3.082',
      'GCTCCTTG': 'R3.083',
      'GTAAGGTG': 'R3.084',
      'GAGGATGG': 'R3.085',
      'GTTGTCGG': 'R3.086',
      'GGATTAGG': 'R3.087',
      'GATAGAGG': 'R3.088',
      'GTGTGTCG': 'R3.089',
      'GCAATCCG': 'R3.090',
      'GACCTTAG': 'R3.091',
      'GCCTGTTC': 'R3.092',
      'GCACTGTC': 'R3.093',
      'GCTAACTC': 'R3.094',
      'GATTCATC': 'R3.095',
      'GTCTTGGC': 'R3.096',
      'ACTGGTAG': 'R3.097',
      'GCTCCAAC': 'R3.098',
      'GCGTAAGA': 'R3.099',
      'TGACCATC': 'R3.100',
      'GGATTCTC': 'R3.101',
      'TCGTTGCG': 'R3.102',
      'GCCAGACT': 'R3.103',
      'TAGGTCCG': 'R3.104',
      'CCGCCTCC': 'R3.105',
      'AGTAGATC': 'R3.106',
      'AGTCTTAA': 'R3.107',
      'GAATGAAT': 'R3.108',
      'CTGAAGTC': 'R3.109',
      'TTATAGCT': 'R3.110',
      'CTTCTCTT': 'R3.111',
      'CCTATTGA': 'R3.112',
      'GCGGCGCA': 'R3.113',
      'AGGCAGTA': 'R3.114',
      'GTTAATTC': 'R3.115',
      'GGAGAATG': 'R3.116',
      'TCTGCCGG': 'R3.117',
      'GCCTCGGT': 'R3.118',
      'CTAATATG': 'R3.119',
      'TGGCTGAC': 'R3.120',
      'ATTGCATA': 'R3.121',
      'AGTCGCGT': 'R3.122',
      'GTTCCGTA': 'R3.123',
      'ATCCGGAT': 'R3.124',
      'TATTACTC': 'R3.125',
      'ACGGACCT': 'R3.126',
      'GTAAGGCC': 'R3.127',
      'GAGACGTT': 'R3.128',
      'TTGGAACT': 'R3.129',
      'AAGGTTAA': 'R3.130',
      'TTAGCGGT': 'R3.131',
      'TATCAAGG': 'R3.132',
      'TGGTATCA': 'R3.133',
      'CTCGATCT': 'R3.134',
      'CTAGCTTC': 'R3.135',
      'CGAGGCCG': 'R3.136',
      'AACCTTCT': 'R3.137',
      'TCGGAGTA': 'R3.138',
      'GCAAGTTA': 'R3.139',
      'GAGCTCAG': 'R3.140',
      'CAGGAATC': 'R3.141',
      'GGTTGGCG': 'R3.142',
      'TGCTCCTT': 'R3.143',
      'CCGGCGTT': 'R3.144',
      'CAAGCATA': 'R3.145',
      'TACGCAGG': 'R3.146',
      'GGCGTCTG': 'R3.147',
      'ATGACGGC': 'R3.148',
      'AGACTACT': 'R3.149',
      'AATCTCGC': 'R3.150',
      'GGCGGATA': 'R3.151',
      'GCGCCTAG': 'R3.152',
      'TGCTTAGC': 'R3.153',
      'TTAGTATA': 'R3.154',
      'TCTAGCAA': 'R3.155',
      'TCCAATGG': 'R3.156',
      'CCAACCGC': 'R3.157',
      'GGTAACTA': 'R3.158',
      'GATGATTG': 'R3.159',
      'ATTATCGA': 'R3.160',
      'TAGATAAT': 'R3.161',
      'CTTACGTT': 'R3.162',
      'TACTATTA': 'R3.163',
      'AGAAGCGA': 'R3.164',
      'ATCTAGTT': 'R3.165',
      'AAGTTCCA': 'R3.166',
      'GAACTAAC': 'R3.167',
      'CCGTCGAA': 'R3.168',
      'TGATTGAA': 'R3.169',
      'TTCTAATC': 'R3.170',
      'GATAGCCG': 'R3.171',
      'TTGACTTC': 'R3.172',
      'GTCGAGGT': 'R3.173',
      'GCTATCCT': 'R3.174',
      'ATAGAGAG': 'R3.175',
      'TTCTTGAC': 'R3.176',
      'GTTAGGAG': 'R3.177',
      'AGTCATGG': 'R3.178',
      'AGACCAGG': 'R3.179',
      'TGCCGGTA': 'R3.180',
      'TAACTGGC': 'R3.181',
      'GTCGGCAA': 'R3.182',
      'CAACGGTC': 'R3.183',
      'ATGCTAAC': 'R3.184',
      'TGGAGGCG': 'R3.185',
      'CCTTGGCT': 'R3.186',
      'GGATGGTA': 'R3.187',
      'GAGGAGAG': 'R3.188',
      'CGGATAAG': 'R3.189',
      'TGCAGAGA': 'R3.190',
      'TTGGCCAG': 'R3.191',
      'TAACCAGT': 'R3.192'}
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


#if i1_in == "NA":
#    print('2 fqs are supplied')
#else:
#    print('4 fqs are supplied')

#if options.qc == True:
#    print('Running QC on ' + str(qcreads+1)  + " reads")
    
# name outputs and print to working dir
p1_file = p1_in.split('/')[-1]
p2_file = p2_in.split('/')[-1]

#check for file type and open input file
append = p1_in.split('.')[-1]
if append == "gz":
    p1_rds = io.BufferedReader(gzip.open(p1_in,'rb'))
    p2_rds = io.BufferedReader(gzip.open(p2_in,'rb'))
    # set up files for undetermined reads
    dis1 = io.BufferedWriter(open(prefix + "discard.R1.fq", 'wb'))
    dis2 = io.BufferedWriter(open(prefix + "discard.R2.fq", 'wb'))
    # read in optional index files
    if i1_in != "NA":
        i1_rds = io.BufferedReader(gzip.open(i1_in,'rb'))
        i2_rds = io.BufferedReader(gzip.open(i2_in,'rb'))
else:
    sys.exit("ERROR! The input file2 must be a .fastq.gz")

##### SCRIPT #####
# initialize variables
i=0;j=0;k=0;tot_b=0;count=1
n=20  # match seq
mismatch=1  # only allow 0-1 mismatches for now
good_r=0

# check if reads are indexed
p1_line = p1_rds.readline()
seqhead1 = p1_line.decode()
# print("first read header:")
# print(seqhead1)

if "+" in seqhead1 or "_" in seqhead1:
    indexed = True
#    print("Fastqs are properly indexed")
else:
    indexed = False
#    print("Will update index")
    if i1_in == "NA":
        sys.exit("warning: one pair of the fastqs are not indexed, but index reads are not supplied")
        

p1_rds.close()
p1_rds = io.BufferedReader(gzip.open(p1_in,'rb'))

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
    if os.path.exists(outName + ".R1.fq.gz"):
        os.remove(outName + ".R1.fq.gz")
    if os.path.exists(outName + ".R2.fq.gz"):
        os.remove(outName + ".R2.fq.gz")

    for primer in metaData['Primer']:
        project[primer] = outName
    for primer in metaData['Primer']:
        sampletype[primer] = outType
#    for primer in metaData['Primer']:
#        N6type[primer] = "F"
#        if sampletype[primer] == 'RNA':
#            if 'N6' in metaData.keys():
#                N6type[primer] = metaData['N6']
            
#print(project)
#print(sampletype)
#print(N6type)

# generate barcode set
r1set = dict()
r2set = dict()
r3set = dict()
p5set = dict()
for barcode, name in R1.items():
    barcodes = barcodeSet(barcode)
    for bc in barcodes:
        r1set[bc] = name
for barcode, name in R2.items():
    barcodes = barcodeSet(barcode)
    for bc in barcodes:
        r2set[bc] = name
for barcode, name in R3.items():
    barcodes = barcodeSet(barcode)
    barcodes2 = barcodeSet(barcode[0:7])
    barcodes3 = barcodeSet(barcode[0:6])
    barcodes4 = barcodeSet(barcode[0:5])
    for bc in barcodes:
        r3set[bc] = name
    for bc in barcodes2:
        r3set[bc] = name
    for bc in barcodes3:
        r3set[bc] = name
    for bc in barcodes4:
        r3set[bc] = name
for barcode, name in P5revcomp.items():
    barcodes = barcodeSet(barcode)
    for bc in barcodes:
        p5set[bc] = name
for barcode, name in P5fwdstr.items():
    barcodes = barcodeSet(barcode)
    for bc in barcodes:
        p5set[bc] = name

files_r1 = dict()
files_r2 = dict()
for proj in project.values():
    f1 = io.BufferedWriter(open((prefix + proj + ".R1.fq"), 'ab'))
    f2 = io.BufferedWriter(open((prefix + proj + ".R2.fq"), 'ab'))
    files_r1[proj] = f1
    files_r2[proj] = f2
#print(files_r1)

start = time.process_time()        
while 1:
    # read lines
    p1_line = p1_rds.readline()
    p2_line = p2_rds.readline()

    # add index to biological reads
    if indexed == False:
        i1_line = i1_rds.readline()
        i2_line = i2_rds.readline()
        
    # break if at end of file
    if not p1_line:
        break

    # load fastq into memory
    if count ==1:
        seqhead1 = p1_line.decode()
        seqhead2 = p2_line.decode()
#        print(seqhead1 + str(i))
        # modify read header
        seqhead1 = seqhead1.replace("1:N:0:1", "1:N:0:")
        seqhead1 = seqhead1.replace("1:N:0:2", "1:N:0:")
        seqhead1 = str.encode(seqhead1.replace("1:N:0:0", "1:N:0:"))
        seqhead2 = seqhead2.replace("4:N:0:1", "2:N:0:")
        seqhead2 = seqhead2.replace("4:N:0:2", "2:N:0:")
        seqhead2 = str.encode(seqhead2.replace("2:N:0:0", "2:N:0:"))        
    elif count ==2:
        seq1 = p1_line.rstrip()
        seq2 = p2_line.rstrip()
        # skip reads in empty
        if seq1.decode() == "+":
            print("warning: found empty biological read. better to use Untrimmed fastq")
            i = i + 1
            p1_line = p1_rds.readline()
            p1_line = p1_rds.readline()
            p2_line = p2_rds.readline()
            p2_line = p2_rds.readline()
            seqhead1 = p1_line.decode()
            seqhead2 = p2_line.decode()
            seqhead1 = seqhead1.replace("1:N:0:1", "1:N:0:")
            seqhead1 = seqhead1.replace("1:N:0:2", "1:N:0:")
            seqhead1 = str.encode(seqhead1.replace("1:N:0:0", "1:N:0:"))
            seqhead2 = seqhead2.replace("4:N:0:1", "2:N:0:")
            seqhead2 = seqhead2.replace("4:N:0:2", "2:N:0:")
            seqhead2 = str.encode(seqhead2.replace("2:N:0:0", "2:N:0:"))
            if indexed == False:
                i1_line = i1_rds.readline()
                i1_line = i1_rds.readline()
                i1_line = i1_rds.readline()
                i1_line = i1_rds.readline()
                i2_line = i2_rds.readline()
                i2_line = i2_rds.readline()
                i2_line = i2_rds.readline()
                i2_line = i2_rds.readline()
#            print(seqhead1)
#            print(i1_line)

        # print(seq1)
        # print(i1_line)
        if indexed == False:
            id1 = i1_line.rstrip()
            id2 = i2_line.rstrip()
            # skip lines in index read is empty
            if "F" in id1.decode():
                print("warning: found empty index read. better to use Untrimmed fastq")
                
            seqhead1 = str.encode(seqhead1.decode().replace("\n", ""))
            seqhead2 = str.encode(seqhead2.decode().replace("\n", ""))
            seqhead1 = (seqhead1 + id1 + b"+" + id2 + b"\n")
            seqhead2 = (seqhead2 + id1 + b"+" + id2 + b"\n")
#        print(id1)    
        # update barcode to R1.xx,R2.xx,R3.xx,P5.xx
        barcode = ""
        barcodeMatch = 0
        index = seqhead1.decode().find("1:N:0:")
        # seqheadMod = seqhead1.decode().replace("_", "+")
        index2 = seqhead1.decode().find("+")
        r1 = seqhead1[(index + 21):(index + 29)].decode()
        if (r1 in r1set):
            barcode = r1set[r1]
            barcodeMatch += 1
        r2 = seqhead1[(index + 59):(index + 67)].decode()
        if (r2 in r2set):
            barcode = barcode + "," + r2set[r2]
            barcodeMatch += 1
        r3 = seqhead1[(index + 97):(index + 105)].decode()
        if (r3 in r3set):
            barcode = barcode + "," + r3set[r3]
            barcodeMatch += 1
        seqheadMod = seqhead1.decode().replace("_", "+")
        index2 = seqheadMod.find("+")
#        p5 = seqhead1[(index2 + 1): -1].decode()
        p5 = seqhead1[(index + 106):(index + 114)].decode()
        if (p5 in p5set):
            barcode = barcode + "," + p5set[p5]
            barcodeMatch += 1
#            print(barcode)
#            print(barcodeMatch)
    elif count ==3:
        qualhead1 = p1_line
        qualhead2 = p2_line
    elif count ==4:
        qual1 = p1_line.rstrip()
        qual2 = p2_line.rstrip()
        needtrim = "F"
        if (barcodeMatch == 4):
            if p5set[p5] in project:
#                print(p5set[p5])
#                print(sampletype[p5set[p5]])
                if sampletype[p5set[p5]] == "ATAC" or sampletype[p5set[p5]] == "TAPS":
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                    needtrim = "T"
                elif sampletype[p5set[p5]] == "DipC":
                    seqhead1 = seqhead1[ :index-1] + b"_1:N:0:" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_1:N:0:" + str.encode(barcode) + b"\n"
                    needtrim = "T"
                elif sampletype[p5set[p5]] == "RNA":
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                elif sampletype[p5set[p5]] == "crop":
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                    index3 = seq1.decode().find("GTTTTAG")
                    # remove 22 bp on 5' and anything after GTTTTAG
                    if index3 == -1:
                        seq1 = seq1[22: ]
                        qual1 = qual1[22: ]
                    else:
                        seq1 = seq1[22: index3]
                        qual1 = qual1[22: index3]
                elif sampletype[p5set[p5]] == "cite":
                    # remove 21 bp on 5' and keep 22-32
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                    seq1 = seq1[21: 31]
                    qual1 = qual1[21: 31]
                elif sampletype[p5set[p5]] == "cellhash":
                    # remove 20 bp on 5' and keep 21-31
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[0:10] + b"\n"
                    seq1 = seq1[20: 30]
                    qual1 = qual1[20: 30]
                elif sampletype[p5set[p5]] == "notrim":
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                else:
                    seqhead1 = seqhead1[ :index] + b"1:N:0:" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index] + b"2:N:0:" + str.encode(barcode) + b"\n"
        # align reads to themselves
        i = i+1  # total reads

        # trim reads
        if needtrim == "T":
            rc_seq2 = reverse_complement(seq2[0:n])
            idx = seq1.rfind(str.encode(rc_seq2)) # look for perfect match
        
            if idx > 0:
                j = j+1  # 0 mismatchs
            elif mismatch>0:
                hold = fuzz_align(rc_seq2,seq1.decode(),mismatch)  # else allow for mismatch
                if hold:
                    idx,mis=hold
                    if mis == 1:
                        k=k+1  # 1 mismatch

            # trim reads if idx exist
            if idx > 0:
                # keep track on how much trimming
                tot_b = tot_b+len(seq2[idx+n:-1]) #track total bases trimmed 
            
                # trim data
                seq1 = seq1[0:idx+n-1]
                # modified to sub1 because some aligners (bowtie) dont like perfectly overlapping reads
                seq2 = seq2[0:idx+n-1]
                qual1 = qual1[0:idx+n-1]
                qual2 = qual2[0:idx+n-1]
        # print data
        if barcodeMatch == 4:
            if p5set[p5] in project:                 
                outName = project[p5set[p5]]
                f1 = files_r1[outName]
                f2 = files_r2[outName]               
                f1.write(seqhead1);f1.write(seq1+b"\n")
                f1.write(qualhead1);f1.write(qual1+b"\n")
                f2.write(seqhead2);f2.write(seq2+b"\n")
                f2.write(qualhead2);f2.write(qual2+b"\n")
            good_r = good_r + 1
        else:
            dis1.write(seqhead1);dis1.write(seq1+b"\n")
            dis1.write(qualhead1);dis1.write(qual1+b"\n")
            dis2.write(seqhead2);dis2.write(seq2+b"\n")
            dis2.write(qualhead2);dis2.write(qual2+b"\n")
        if options.qc == True:
            if i > qcreads:
                break
    # increment count
    count = count + 1
    if count == 5:
        count = 1
    else:
        count = count

# close files to write the file
for f in files_r1.values():
    f.close()
for f in files_r2.values():
    f.close()
p1_rds.close();p2_rds.close()
dis1.close();dis2.close()
time = (time.process_time() - start)/60

# print("%.2g" % time + " min comsumed")

## give summary
try:
#    print(str(i)+" sequences total")
#    print(str(j)+" sequences trimmed with 0 mismatches")
#    print(str(k)+" sequences trimmed with 1 mismatch")
#    mean = tot_b/(j+k)
#    print("%.3g" % mean +" mean number of bases trimmed for reads requiring trimming")
    print(str(good_r)+" sucessfully demultiplexed for "+prefix)
    perc = good_r/(i)*100
    print("%.3g" % perc + "% reads demultiplexed")
except ZeroDivisionError: 
    print("Warning too few reads")
