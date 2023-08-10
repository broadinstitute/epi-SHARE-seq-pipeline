#!/usr/bin/env python3

# This script is used to infer molecular barcodes
# from raw sequencing BCL data.
#
# It requires running Picard ExtractIlluminaBarcodes with BARCODE=N,
# to extract all barcodes into *_barcode.txt.gz files first.

import glob
import gzip
import sys

[
    _name,
    multiplex_params_file,
    candidate_molecular_barcodes_file,
    barcode_matches_file,
] = sys.argv

YIELD_THRESHOLD = 0.1
MIN_READ_COUNT = 1e6


def parse_barcodes(file_path):
    with open(file_path) as f:
        barcodes = {}
        for row in f.readlines():
            row = row.strip().split('\t')
            copa = row[0]
            barcode = ''.join(row[1:])
            barcodes[barcode] = copa
        return barcodes


copa_barcodes = parse_barcodes(multiplex_params_file)
molecular_barcodes = parse_barcodes(candidate_molecular_barcodes_file)

# count each unique barcode combination
counts = {}
for extracted in glob.glob('*_barcode.txt.gz'):
    with gzip.open(extracted, 'rt') as f:
        for row in f.readlines():
            barcode = row.split('\t')[0]
            if barcode in counts:
                counts[barcode] += 1
            else:
                counts[barcode] = 1

# add any missing barcodes from the list of CoPAs
for barcode in copa_barcodes.keys():
    if barcode not in counts:
        counts[barcode] = 0


def distance(b1, b2):
    return sum(c1 != c2 for c1, c2 in zip(b1, b2))


COPA_UNDEFINED = 'UNDEFINED'

# match barcodes to candidates
results = {}
molecular_barcode_len = len(next(iter(molecular_barcodes)))
for barcode, count in counts.items():
    molecular_barcode_matched = False
    molecular_barcode_match = barcode[:molecular_barcode_len]
    molecular_barcode_match_name = molecular_barcode_match

    if molecular_barcode_match in molecular_barcodes:
        molecular_barcode_matched = True
        molecular_barcode_match_name = molecular_barcodes[molecular_barcode_match]

    barcode_match = molecular_barcode_match
    copa = copa_barcodes[barcode_match] if barcode_match in copa_barcodes else COPA_UNDEFINED
    if barcode_match in results:
        results[barcode_match]['Count'] += count
    else:
        results[barcode_match] = {
            'CoPA': copa,
            'Molecular Barcode': molecular_barcode_match_name,
            'Count': count,
            'Matched': molecular_barcode_matched
        }

# show barcodes that correspond to a CoPA or have a matched
# barcode at the top of the output file, otherwise
# sort by count
results = sorted(
    results.values(),
    key=lambda r: (r['CoPA'], int(not r['Matched']), -r['Count'])
)

# calculate % of average yield
total_yield = 0
copa_count = 0
for r in results:
    if r['CoPA'] != COPA_UNDEFINED:
        total_yield += r['Count']
        copa_count += 1
avg_yield = total_yield / copa_count if copa_count else None
for r in results:
    percent_avg_yield = ''
    if r['CoPA'] != COPA_UNDEFINED:
        percent_avg_yield = '{:.2f}%'.format(
            100 * r['Count'] / avg_yield) if avg_yield else 0
    r['Percent of average'] = percent_avg_yield

# report results as a TSV
with open(barcode_matches_file, 'w') as f:
    header = (
        'CoPA', 'Molecular Barcode',
        'Count', 'Percent of average',
    )
    print('\t'.join(header), file=f)

    # print CoPA matches, barcode matches, and the top barcodes
    # with the highest read count
    for r in results:
        if (
            r['CoPA'] != COPA_UNDEFINED or
            r['Matched'] or
            r['Count'] >= MIN_READ_COUNT
        ):
            print('\t'.join((str(r[col]) for col in header)), file=f)

# fail the task (and the workflow) for low yield
if not avg_yield:
    raise Exception('None of the candidate barcodes matched any CoPAs!')
failed_copas = []
for r in results:
    if (r['CoPA'] != COPA_UNDEFINED and
            float(r['Percent of average'].replace('%', '')) < YIELD_THRESHOLD):
        failed_copas.append(r['CoPA'])
failed_copas = ', '.join(failed_copas)
if failed_copas:
    raise Exception(
        f'Found CoPA(s) with < {YIELD_THRESHOLD}% yield: {failed_copas}')