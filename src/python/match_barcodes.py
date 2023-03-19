import gzip

import numpy as np
# import pandas as pd

import matcha
import sys

REV_COMP = str.maketrans("ATGC", "TACG")
def reverse_complement(seq):
    return str.translate(seq, REV_COMP)[::-1]

def get_open_fn(path):
    with open(path, "rb") as f:
        is_gzipped = (f.read(2) == b'\x1f\x8b')
    return gzip.open if is_gzipped else open

def read_barcodes(path, revcomp):
    # if path.endswith(".tsv"):
    #     bc = pd.read_csv(path, sep="\t")["sequence"]
    # else:
    open_fn = get_open_fn(path)
    with open_fn(path, 'rt') as file:
        bc = [b.strip() for b in file]
    if revcomp:
        valid = [reverse_complement(b) for b in bc]
    else:
        valid = bc

    return valid

def match_one_bc(fastqs, whitelists, revcomp, max_barcode_dist, offsets, fastq1_out_path, fastq2_out_path, qc_path, threads):
    f = matcha.FastqReader(threads = threads)
    f.add_sequence("R1", fastqs["R1"], output_path=fastq1_out_path)
    f.add_sequence("R2", fastqs["R2"])
    f.add_sequence("R3", fastqs["R3"], output_path=fastq2_out_path)

    with open(revcomp["R2"]) as rf:
        rc = (int(rf.read().strip()) == 1)

    barcode_sequences = read_barcodes(whitelists["R2"], rc)
    cell_barcode = matcha.HashMatcher(
        sequences = barcode_sequences,
        labels = barcode_sequences,
        max_mismatches=max_barcode_dist,
        subsequence_count=2
    )
    f.add_barcode("cell", cell_barcode, "R2", match_start=offsets["R2"])
    f.set_output_names("{read_name} CB:Z:{cell}")

    barcode_counts = np.zeros(max_barcode_dist + 2, int)

    total_reads = 0
    total_pass = 0

    # print("start read") ####
    chunk_size = 10000
    while f.read_chunk(chunk_size):
        pass_filter = (f.get_match_result("cell", "dist") <=max_barcode_dist) & \
            (f.get_match_result("cell", "second_best_dist") > f.get_match_result("cell", "dist"))

        total_reads += len(pass_filter)
        total_pass += pass_filter.sum()
        values, counts = np.unique(f.get_match_result("cell", "dist"), return_counts=True)
        barcode_counts[np.minimum(values, max_barcode_dist + 1)] += counts

        f.write_chunk(pass_filter)

    with open(qc_path, "w") as stats_output:
        print(f"{total_pass}/{total_reads} reads passing, ({total_pass/total_reads*100:.2f}%)\n", file=stats_output)
        print("mismatches\treads", file=stats_output)
        for dist in range(max_barcode_dist + 2):
            print(
                dist if dist <= max_barcode_dist else f">{max_barcode_dist}",
                barcode_counts[dist],
                sep = "\t",
                file=stats_output
            )


# def match_two_bc(fastqs, whitelists, revcomp, max_barcode_dist, offsets, fastq1_out_path, fastq2_out_path, qc_path, threads):
#     f = matcha.FastqReader(threads = threads)
#     f.add_sequence("R1", fastqs["R1"], output_path=fastq1_out_path)
#     f.add_sequence("R2", fastqs["R2"], output_path=fastq2_out_path)
#     f.add_sequence("I1", fastqs["I1"])
#     f.add_sequence("I2", fastqs["I2"])

#     i5_sequences, i5_maybe_rc = read_barcodes(whitelists["I2"], revcomp["I2"])
#     T7_sequences, T7_maybe_rc = read_barcodes(whitelists["I1"], revcomp["I1"])

#     i5_barcode = matcha.HashMatcher(
#         sequences = i5_maybe_rc,
#         labels = i5_sequences,
#         max_mismatches=max_barcode_dist,
#         subsequence_count=2
#     )

#     T7_barcode = matcha.HashMatcher(
#         sequences = T7_maybe_rc,
#         labels = T7_sequences,
#         max_mismatches=max_barcode_dist,
#         subsequence_count=2
#     )

#     f.add_barcode("i5", i5_barcode, "I2", match_start=offsets["I2"])
#     f.add_barcode("T7", T7_barcode, "I1", match_start=offsets["I1"])

#     f.set_output_names("{read_name} CB:Z:{i5}{T7}")

#     barcode_counts = np.zeros((max_barcode_dist + 2, max_barcode_dist + 2), int)

#     total_reads = 0
#     total_pass = 0

#     chunk_size = 10000

#     dists = [None, None]
#     second_dists = [None, None]
#     while f.read_chunk(chunk_size):
#         dists[0] = f.get_match_result("i5", "dist")
#         second_dists[0] = f.get_match_result("i5", "second_best_dist")
#         dists[1] = f.get_match_result("T7", "dist")
#         second_dists[1] = f.get_match_result("T7", "second_best_dist")

#         pass_filter = (dists[0] < max_barcode_dist) & \
#             (dists[1] < max_barcode_dist) & \
#             (dists[0] + dists[1] < second_dists[0] + second_dists[1])

#         total_reads += len(pass_filter)
#         total_pass += pass_filter.sum()

#         values, counts = np.unique(dists, axis = 1, return_counts=True)
#         indices = np.minimum(values, max_barcode_dist+1)
#         barcode_counts[(indices[0], indices[1])] += counts

#         f.write_chunk(pass_filter)

#     with open(qc_path, "w") as stats_output:
#         print(f"{total_pass}/{total_reads} reads passing, ({total_pass/total_reads*100:.2f}%)\n", file=stats_output)
#         print("mismatches_i5\tmismatches_T7\treads", file=stats_output)
#         for i5_dist in range(max_barcode_dist + 2):
#             for T7_dist in range(max_barcode_dist + 2):
#                 print(
#                     i5_dist if i5_dist <= max_barcode_dist else f">{max_barcode_dist}",
#                     T7_dist if T7_dist <= max_barcode_dist else f">{max_barcode_dist}",
#                     barcode_counts[i5_dist, T7_dist],
#                     sep = "\t",
#                     file=stats_output
#                 )

modality = sys.argv[4]
whitelist = sys.argv[7]
fastq1_out_path = sys.argv[8]
fastq2_out_path = sys.argv[9]
qc_path = sys.argv[10]
threads = int(sys.argv[11])
max_barcode_dist = int(sys.argv[5])
fastqs = {
    "R1": sys.argv[1],
    "R2": sys.argv[3],
    "R3": sys.argv[2],
}
revcomp = {
    "R2": sys.argv[6],
}
if modality == "10x":
    whitelists = {
        "R2": whitelist,
    }
    offsets = {
        "R2": 0,
    }
    match_one_bc(fastqs, whitelists, revcomp, max_barcode_dist, offsets, fastq1_out_path, fastq2_out_path, qc_path, threads)

elif modality == "10x_multiome":
    whitelists = {
        "R2": whitelist,
    }
    offsets = {
        "R2": 8,
    }
    match_one_bc(fastqs, whitelists, revcomp, max_barcode_dist, offsets, fastq1_out_path, fastq2_out_path, qc_path, threads)
