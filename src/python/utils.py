#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions
"""


def check_putative_barcode(barcode_str, barcode_exact_dict, barcode_mismatch_dict):
    """
    Check exact match of barcode, then 1 mismatch. Return corrected
    barcode and boolean indicating whether the correction was an exact match.
    If no match was found, return None for the corrected barcode.
    """
    corrected = barcode_exact_dict.get(barcode_str)  # check exact location first

    if corrected:
        return corrected, True
    
    corrected = barcode_mismatch_dict.get(barcode_str)  # check mismatch

    if corrected:
        return corrected, False
    
    return None, False

def create_barcode_dicts(barcode_list):
    """
    Creates dictionaries containing exact match and 1-mismatch sequences
    """
    barcode_exact_dict = dict()  # [seq: seq]
    barcode_mismatch_dict = dict()  # [mismatch: seq]
    for barcode in barcode_list:
        barcode_exact_dict[barcode] = barcode  # exact match
        for i, base in enumerate(barcode):
            for x in 'ACGT':
                if base != x:
                    # add mismatch possibilities at pos i
                    mismatch = barcode[:i] + x + barcode[i + 1:]
                    barcode_mismatch_dict[mismatch] = barcode
    return barcode_exact_dict, barcode_mismatch_dict
