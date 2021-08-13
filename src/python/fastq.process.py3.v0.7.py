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
                    # allow 1 mismatch or 1 bp shift
                    barcodeSet.add((barcode[1:] + base))
                    barcodeSet.add((base + barcode[:-1]))
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
R1 = {'ATCACGTT': 'R1.01',
      'CGATGTTT': 'R1.02',
      'TTAGGCAT': 'R1.03',
      'TGACCACT': 'R1.04',
      'ACAGTGGT': 'R1.05',
      'GCCAATGT': 'R1.06',
      'CAGATCTG': 'R1.07',
      'ACTTGATG': 'R1.08',
      'GATCAGCG': 'R1.09',
      'TAGCTTGT': 'R1.10',
      'GGCTACAG': 'R1.11',
      'CTTGTACT': 'R1.12',
      'TGGTTGTT': 'R1.13',
      'TCTCGGTT': 'R1.14',
      'TAAGCGTT': 'R1.15',
      'TCCGTCTT': 'R1.16',
      'TGTACCTT': 'R1.17',
      'TTCTGTGT': 'R1.18',
      'TCTGCTGT': 'R1.19',
      'TTGGAGGT': 'R1.20',
      'TCGAGCGT': 'R1.21',
      'TGATACGT': 'R1.22',
      'TGCATAGT': 'R1.23',
      'TTGACTCT': 'R1.24',
      'TGCGATCT': 'R1.25',
      'TTCCTGCT': 'R1.26',
      'TAGTGACT': 'R1.27',
      'TACAGGAT': 'R1.28',
      'TCCTCAAT': 'R1.29',
      'TGTGGTTG': 'R1.30',
      'TACTAGTC': 'R1.31',
      'TTCCATTG': 'R1.32',
      'TCGAAGTG': 'R1.33',
      'TAACGCTG': 'R1.34',
      'TTGGTATG': 'R1.35',
      'TGAACTGG': 'R1.36',
      'TACTTCGG': 'R1.37',
      'TCTCACGG': 'R1.38',
      'TCAGGAGG': 'R1.39',
      'TAAGTTCG': 'R1.40',
      'TCCAGTCG': 'R1.41',
      'TGTATGCG': 'R1.42',
      'TCATTGAG': 'R1.43',
      'TGGCTCAG': 'R1.44',
      'TATGCCAG': 'R1.45',
      'TCAGATTC': 'R1.46',
      'TAGTCTTG': 'R1.47',
      'TTCAGCTC': 'R1.48',
      'TGTCTATC': 'R1.49',
      'TATGTGGC': 'R1.50',
      'TTACTCGC': 'R1.51',
      'TCGTTAGC': 'R1.52',
      'TACCGAGC': 'R1.53',
      'TGTTCTCC': 'R1.54',
      'TTCGCACC': 'R1.55',
      'TTGCGTAC': 'R1.56',
      'TCTACGAC': 'R1.57',
      'TGACAGAC': 'R1.58',
      'TAGAACAC': 'R1.59',
      'TCATCCTA': 'R1.60',
      'TGCTGATA': 'R1.61',
      'TAGACGGA': 'R1.62',
      'TGTGAAGA': 'R1.63',
      'TCTCTTCA': 'R1.64',
      'TTGTTCCA': 'R1.65',
      'TGAAGCCA': 'R1.66',
      'TACCACCA': 'R1.67',
      'TGCGTGAA': 'R1.68',
      'GGTGAGTT': 'R1.69',
      'GATCTCTT': 'R1.70',
      'GTGTCCTT': 'R1.71',
      'GACGGATT': 'R1.72',
      'GCAACATT': 'R1.73',
      'GGTCGTGT': 'R1.74',
      'GAATCTGT': 'R1.75',
      'GTACATCT': 'R1.76',
      'GAGGTGCT': 'R1.77',
      'GCATGGCT': 'R1.78',
      'GTTAGCCT': 'R1.79',
      'GTCGCTAT': 'R1.80',
      'GGAATGAT': 'R1.81',
      'GAGCCAAT': 'R1.82',
      'GCTCCTTG': 'R1.83',
      'GTAAGGTG': 'R1.84',
      'GAGGATGG': 'R1.85',
      'GTTGTCGG': 'R1.86',
      'GGATTAGG': 'R1.87',
      'GATAGAGG': 'R1.88',
      'GTGTGTCG': 'R1.89',
      'GCAATCCG': 'R1.90',
      'GACCTTAG': 'R1.91',
      'GCCTGTTC': 'R1.92',
      'GCACTGTC': 'R1.93',
      'GCTAACTC': 'R1.94',
      'GATTCATC': 'R1.95',
      'GTCTTGGC': 'R1.96'}
R2 = {'ATCACGTT': 'R2.01',
      'CGATGTTT': 'R2.02',
      'TTAGGCAT': 'R2.03',
      'TGACCACT': 'R2.04',
      'ACAGTGGT': 'R2.05',
      'GCCAATGT': 'R2.06',
      'CAGATCTG': 'R2.07',
      'ACTTGATG': 'R2.08',
      'GATCAGCG': 'R2.09',
      'TAGCTTGT': 'R2.10',
      'GGCTACAG': 'R2.11',
      'CTTGTACT': 'R2.12',
      'TGGTTGTT': 'R2.13',
      'TCTCGGTT': 'R2.14',
      'TAAGCGTT': 'R2.15',
      'TCCGTCTT': 'R2.16',
      'TGTACCTT': 'R2.17',
      'TTCTGTGT': 'R2.18',
      'TCTGCTGT': 'R2.19',
      'TTGGAGGT': 'R2.20',
      'TCGAGCGT': 'R2.21',
      'TGATACGT': 'R2.22',
      'TGCATAGT': 'R2.23',
      'TTGACTCT': 'R2.24',
      'TGCGATCT': 'R2.25',
      'TTCCTGCT': 'R2.26',
      'TAGTGACT': 'R2.27',
      'TACAGGAT': 'R2.28',
      'TCCTCAAT': 'R2.29',
      'TGTGGTTG': 'R2.30',
      'TACTAGTC': 'R2.31',
      'TTCCATTG': 'R2.32',
      'TCGAAGTG': 'R2.33',
      'TAACGCTG': 'R2.34',
      'TTGGTATG': 'R2.35',
      'TGAACTGG': 'R2.36',
      'TACTTCGG': 'R2.37',
      'TCTCACGG': 'R2.38',
      'TCAGGAGG': 'R2.39',
      'TAAGTTCG': 'R2.40',
      'TCCAGTCG': 'R2.41',
      'TGTATGCG': 'R2.42',
      'TCATTGAG': 'R2.43',
      'TGGCTCAG': 'R2.44',
      'TATGCCAG': 'R2.45',
      'TCAGATTC': 'R2.46',
      'TAGTCTTG': 'R2.47',
      'TTCAGCTC': 'R2.48',
      'TGTCTATC': 'R2.49',
      'TATGTGGC': 'R2.50',
      'TTACTCGC': 'R2.51',
      'TCGTTAGC': 'R2.52',
      'TACCGAGC': 'R2.53',
      'TGTTCTCC': 'R2.54',
      'TTCGCACC': 'R2.55',
      'TTGCGTAC': 'R2.56',
      'TCTACGAC': 'R2.57',
      'TGACAGAC': 'R2.58',
      'TAGAACAC': 'R2.59',
      'TCATCCTA': 'R2.60',
      'TGCTGATA': 'R2.61',
      'TAGACGGA': 'R2.62',
      'TGTGAAGA': 'R2.63',
      'TCTCTTCA': 'R2.64',
      'TTGTTCCA': 'R2.65',
      'TGAAGCCA': 'R2.66',
      'TACCACCA': 'R2.67',
      'TGCGTGAA': 'R2.68',
      'GGTGAGTT': 'R2.69',
      'GATCTCTT': 'R2.70',
      'GTGTCCTT': 'R2.71',
      'GACGGATT': 'R2.72',
      'GCAACATT': 'R2.73',
      'GGTCGTGT': 'R2.74',
      'GAATCTGT': 'R2.75',
      'GTACATCT': 'R2.76',
      'GAGGTGCT': 'R2.77',
      'GCATGGCT': 'R2.78',
      'GTTAGCCT': 'R2.79',
      'GTCGCTAT': 'R2.80',
      'GGAATGAT': 'R2.81',
      'GAGCCAAT': 'R2.82',
      'GCTCCTTG': 'R2.83',
      'GTAAGGTG': 'R2.84',
      'GAGGATGG': 'R2.85',
      'GTTGTCGG': 'R2.86',
      'GGATTAGG': 'R2.87',
      'GATAGAGG': 'R2.88',
      'GTGTGTCG': 'R2.89',
      'GCAATCCG': 'R2.90',
      'GACCTTAG': 'R2.91',
      'GCCTGTTC': 'R2.92',
      'GCACTGTC': 'R2.93',
      'GCTAACTC': 'R2.94',
      'GATTCATC': 'R2.95',
      'GTCTTGGC': 'R2.96'}
R3 = {'ATCACGTT': 'R3.01',
      'CGATGTTT': 'R3.02',
      'TTAGGCAT': 'R3.03',
      'TGACCACT': 'R3.04',
      'ACAGTGGT': 'R3.05',
      'GCCAATGT': 'R3.06',
      'CAGATCTG': 'R3.07',
      'ACTTGATG': 'R3.08',
      'GATCAGCG': 'R3.09',
      'TAGCTTGT': 'R3.10',
      'GGCTACAG': 'R3.11',
      'CTTGTACT': 'R3.12',
      'TGGTTGTT': 'R3.13',
      'TCTCGGTT': 'R3.14',
      'TAAGCGTT': 'R3.15',
      'TCCGTCTT': 'R3.16',
      'TGTACCTT': 'R3.17',
      'TTCTGTGT': 'R3.18',
      'TCTGCTGT': 'R3.19',
      'TTGGAGGT': 'R3.20',
      'TCGAGCGT': 'R3.21',
      'TGATACGT': 'R3.22',
      'TGCATAGT': 'R3.23',
      'TTGACTCT': 'R3.24',
      'TGCGATCT': 'R3.25',
      'TTCCTGCT': 'R3.26',
      'TAGTGACT': 'R3.27',
      'TACAGGAT': 'R3.28',
      'TCCTCAAT': 'R3.29',
      'TGTGGTTG': 'R3.30',
      'TACTAGTC': 'R3.31',
      'TTCCATTG': 'R3.32',
      'TCGAAGTG': 'R3.33',
      'TAACGCTG': 'R3.34',
      'TTGGTATG': 'R3.35',
      'TGAACTGG': 'R3.36',
      'TACTTCGG': 'R3.37',
      'TCTCACGG': 'R3.38',
      'TCAGGAGG': 'R3.39',
      'TAAGTTCG': 'R3.40',
      'TCCAGTCG': 'R3.41',
      'TGTATGCG': 'R3.42',
      'TCATTGAG': 'R3.43',
      'TGGCTCAG': 'R3.44',
      'TATGCCAG': 'R3.45',
      'TCAGATTC': 'R3.46',
      'TAGTCTTG': 'R3.47',
      'TTCAGCTC': 'R3.48',
      'TGTCTATC': 'R3.49',
      'TATGTGGC': 'R3.50',
      'TTACTCGC': 'R3.51',
      'TCGTTAGC': 'R3.52',
      'TACCGAGC': 'R3.53',
      'TGTTCTCC': 'R3.54',
      'TTCGCACC': 'R3.55',
      'TTGCGTAC': 'R3.56',
      'TCTACGAC': 'R3.57',
      'TGACAGAC': 'R3.58',
      'TAGAACAC': 'R3.59',
      'TCATCCTA': 'R3.60',
      'TGCTGATA': 'R3.61',
      'TAGACGGA': 'R3.62',
      'TGTGAAGA': 'R3.63',
      'TCTCTTCA': 'R3.64',
      'TTGTTCCA': 'R3.65',
      'TGAAGCCA': 'R3.66',
      'TACCACCA': 'R3.67',
      'TGCGTGAA': 'R3.68',
      'GGTGAGTT': 'R3.69',
      'GATCTCTT': 'R3.70',
      'GTGTCCTT': 'R3.71',
      'GACGGATT': 'R3.72',
      'GCAACATT': 'R3.73',
      'GGTCGTGT': 'R3.74',
      'GAATCTGT': 'R3.75',
      'GTACATCT': 'R3.76',
      'GAGGTGCT': 'R3.77',
      'GCATGGCT': 'R3.78',
      'GTTAGCCT': 'R3.79',
      'GTCGCTAT': 'R3.80',
      'GGAATGAT': 'R3.81',
      'GAGCCAAT': 'R3.82',
      'GCTCCTTG': 'R3.83',
      'GTAAGGTG': 'R3.84',
      'GAGGATGG': 'R3.85',
      'GTTGTCGG': 'R3.86',
      'GGATTAGG': 'R3.87',
      'GATAGAGG': 'R3.88',
      'GTGTGTCG': 'R3.89',
      'GCAATCCG': 'R3.90',
      'GACCTTAG': 'R3.91',
      'GCCTGTTC': 'R3.92',
      'GCACTGTC': 'R3.93',
      'GCTAACTC': 'R3.94',
      'GATTCATC': 'R3.95',
      'GTCTTGGC': 'R3.96'}
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
    dis1 = io.BufferedWriter(open("discard/discard.R1.fq", 'wb'))
    dis2 = io.BufferedWriter(open("discard/discard.R2.fq", 'wb'))
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
    if os.path.exists("out/" + outName + ".R1.fq.gz"):
        os.remove("out/" + outName + ".R1.fq.gz")
    if os.path.exists("out/" + outName + ".R2.fq.gz"):
        os.remove("out/" + outName + ".R2.fq.gz")

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

file_names = open("path_to_files.txt","w")
file_name_arr1 = []
file_name_arr2 = []

files_r1 = dict()
files_r2 = dict()
for proj in project.values():
    f1 = io.BufferedWriter(open(("out/" + proj + ".R1.fq"), 'ab'))
    f2 = io.BufferedWriter(open(("out/" + proj + ".R2.fq"), 'ab'))
    files_r1[proj] = f1
    files_r2[proj] = f2
    file_name_arr1.append("out/" + proj + ".R1.fq")
    file_name_arr2.append("out/" + proj + ".R2.fq")

for file1,file2 in set(zip(file_name_arr1,file_name_arr2)):
    file_names.write(f"{file1}\n")
    file_names.write(f"{file2}\n")
#print(files_r1)
file_names.close()


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
                    seqhead1 = seqhead1[ :index-1] + b"_" + seq2[0:10] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                elif sampletype[p5set[p5]] == "crop":
                    seqhead1 = seqhead1[ :index-1] + b"_" + seq2[0:10] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
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
                    seqhead1 = seqhead1[ :index-1] + b"_" + seq2[0:10] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                    seq1 = seq1[21: 31]
                    qual1 = qual1[21: 31]
                elif sampletype[p5set[p5]] == "cellhash":
                    # remove 20 bp on 5' and keep 21-31
                    seqhead1 = seqhead1[ :index-1] + b"_" + seq2[0:10] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
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

        # trimm reads
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
    print(str(good_r)+" sucessfully demultiplexed")
    perc = good_r/(i)*100
    print("%.3g" % perc + "% reads demultiplexed")
except ZeroDivisionError: 
    print("Warning too few reads")
