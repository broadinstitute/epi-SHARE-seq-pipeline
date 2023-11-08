#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions
"""


def check_putative_barcode(barcode_str, quality_str, barcode_exact_dict, barcode_mismatch_dict):
    """
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right
    shift
    """
    correction_type = None
    corrected = barcode_exact_dict.get(barcode_str[1:9])  # check exact location first
    quality = quality_str[1:9]
    if corrected:
        correction_type = "E"
    else:
        corrected = barcode_mismatch_dict.get(barcode_str[1:9])  # check mismatch
        if corrected:
            correction_type = "M"
        else:
            corrected = barcode_mismatch_dict.get(barcode_str[:8])  # check 1bp shift left
            quality = quality_str[:8]
            if corrected:
                correction_type = "L"
            else:  # check 1bp shift right; round 3 is shorter so add "N" for those
                if len(barcode_str) < 10:
                    corrected = barcode_exact_dict.get(barcode_str[2:]+"N")
                    quality = quality_str[2:]+"F"
                    if corrected:
                        correction_type = "R"
                else:
                    corrected = barcode_exact_dict.get(barcode_str[2:])
                    quality = quality_str[2:]
                    if corrected:
                        correction_type = "R"
    return corrected, quality, correction_type

def create_barcode_dicts(barcode_list):
    """
    Creates dictionaries containing exact match and 1-mismatch sequences
    """
    barcode_exact_dict = dict()  # [seq: seq]
    barcode_mismatch_dict = dict()  # [mismatch: seq]
    for barcode in barcode_list:
        barcode_exact_dict[barcode] = barcode  # exact match
        for i, base in enumerate(barcode):
            for x in 'ACGTN':
                if base != x:
                    # add mismatch possibilities at pos i
                    mismatch = barcode[:i] + x + barcode[i + 1:]
                    barcode_mismatch_dict[mismatch] = barcode
    return barcode_exact_dict, barcode_mismatch_dict
