import time
import itertools
import gzip
import heapq
import sys
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
# from kneed import KneeLocator

"""
From https://github.com/kundajelab/ENCODE_scatac/blob/master/workflow/scripts/detect_multiplets.py
"""



def print_and_log(text, outfile, starttime=0):
    logtime = time.process_time() - starttime
    if logtime < 60:
        logtime = "{:,}s".format(logtime)
    else:
        logtime = "{:,}m {:,}s".format(logtime // 60, logtime % 60)
    outfile.write("{} - {}\n".format(logtime, text))
    print("{} - {}".format(logtime, text))

def multiplet_fdr(samples, nulls, fdr_thresh, min_cutoff):
    null_total = nulls.shape[0]
    sample_total = samples.shape[0]

    p = 1 - np.searchsorted(nulls, samples) / null_total
    # p_corr = 1 - (1 - p)**max_beads_per_drop
    q = (p * sample_total) / (sample_total - np.arange(sample_total))
    candidiates, = np.nonzero(q <= fdr_thresh)
    if candidiates.size == 0:
        cut = min_cutoff
    else:
        cut = max(min_cutoff, samples[candidiates.min()])

    return cut, q

def plot_dist(cut, q, samples, nulls, title, x_label, out_path, log_x=False, hist_bins=200):
    fig, ax = plt.subplots(tight_layout=True)
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    if log_x:
        ax.set_xscale('log')
        hist_bins = np.geomspace(min(samples[0], nulls[0]), max(samples[-1], nulls[-1]), hist_bins)

    ax.hist(nulls, bins=hist_bins, alpha=0.5, color="k")
    ax.hist(samples, bins=hist_bins, alpha=0.5, color="b")
    ax2.plot(samples, q, color="g")
    ax.axvline(x=cut, color="r")

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Histogram Frequency")
    ax2.set_ylabel("Non-Monotonicized Q-Value")

    plt.savefig(out_path)

def main(fragments, barcodes_strict, barcodes_expanded, summary, barcodes_status, jac_plot, min_counts=500, fdr_thresh=0.1, max_beads_per_drop=6, min_cutoff=1e-3):
    logout = open(summary, "w")
    starttime = time.process_time()

    print_and_log("Identifying candidate barcodes", logout, starttime)

    cur_clique = set()
    cur_coord = None
    i = 0
    barcode_counts = Counter()
    with gzip.open(fragments, 'rt') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            chr, start, end, barcode = line[:4]

            this_coord = (chr, start, end)
            if this_coord != cur_coord:
                if len(cur_clique) <= max_beads_per_drop:
                    for b in cur_clique:
                        barcode_counts[b] += 1

                cur_clique = set([barcode])
                cur_coord = this_coord

            else:
                cur_clique.add(barcode)

            i += 1
            if i%1e7==0:
                print(i)


    print_and_log(
        f"Original run had {len(barcode_counts)} total cell barcodes",
        logout,
        starttime,
    )

    barcodes_considered = set(k for k, v in barcode_counts.items() if v >= min_counts)
    num_bc = len(barcodes_considered)

    print_and_log(
        f"Identified {num_bc} total barcodes for multiplet detection",
        logout,
        starttime,
    )

    print_and_log("Reading fragments", logout, starttime)

    pair_counts = Counter()
    cur_clique = set()
    cur_coord = None
    i = 0
    with gzip.open(fragments, 'rt') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            chr, start, end, barcode = line[:4]

            if barcode not in barcodes_considered:
                continue

            this_coord = (chr, start, end)
            if this_coord != cur_coord:
                if len(cur_clique) <= max_beads_per_drop:
                    for x, y in itertools.combinations(cur_clique, 2):
                        x, y = (x, y) if x < y else (y, x)
                        pair_counts[(x, y)] += 1

                cur_clique = set([barcode])
                cur_coord = this_coord

            else:
                cur_clique.add(barcode)

            i += 1
            if i%1e7==0:
                print(i)

    print_and_log("Identifying barcode multiplets", logout, starttime)

    print_and_log(
        f"Considered {len(pair_counts)} barcode pairs",
        logout,
        starttime,
    )

    expanded_data = {}
    jac_dists_max = {}
    jac_dists_top = {}
    top_len = max_beads_per_drop + 1
    for x, y in pair_counts.items():
        a, b = x
        bca = barcode_counts[a]
        bcb = barcode_counts[b]
        jac = y/(bca + bcb - y)
        data = [a, b, bca, bcb, y, jac, None]
        expanded_data[x] = data
        if jac > 0:
            jac_dists_max[a] = max(jac_dists_max.get(a, 0), jac)
            jac_dists_max[b] = max(jac_dists_max.get(b, 0), jac)

            aheap = jac_dists_top.setdefault(a, [])
            if len(aheap) < top_len:
                heapq.heappush(aheap, jac)
            else:
                heapq.heappushpop(aheap, jac)

            bheap = jac_dists_top.setdefault(b, [])
            if len(bheap) < top_len:
                heapq.heappush(bheap, jac)
            else:
                heapq.heappushpop(bheap, jac)

    with gzip.open(barcodes_expanded, 'wt') as f:
        f.write("Barcode1\tBarcode2\tBarcode1Counts\tBarcode2Counts\tCommon\tJaccardIndex\n")
        for x, data in expanded_data.items():
            a, b = x
            f.write("{}\t{}\t{}\t{}\t{}\t{:.4f}\n".format(*data[:-1]))

    jac_dists_ref_filt = {}
    jac_dists_max_filt = {}
    jac_dists_ratios = []
    for k, v in jac_dists_top.items():
        if len(v) < top_len:
            continue
        dist_max = jac_dists_max[k]
        dist_ref = v[0]
        # if dist_ref < (1 / min_counts):
        #     continue
        jac_dists_ref_filt[k] = dist_ref
        jac_dists_max_filt[k] = dist_max
        jac_dists_ratios.append(dist_max/dist_ref)

    jac_shift = np.median(jac_dists_ratios)

    samples = np.fromiter(jac_dists_max_filt.values(), dtype=float, count=len(jac_dists_max_filt))
    samples.sort()
    nulls = np.fromiter(jac_dists_ref_filt.values(), dtype=float, count=len(jac_dists_ref_filt)) * jac_shift
    nulls.sort()

    cut, q = multiplet_fdr(samples, nulls, fdr_thresh, min_cutoff)
    plot_dist(cut, q, samples, nulls, "Multiplet Thresholding", "Max Marginal Jaccard Distance", jac_plot, log_x=True)

    min_jac = cut

    print_and_log(
        f"Setting multiplet threshold as {min_jac} for minimum pairwise Jaccard distance",
        logout,
        starttime,
    )

    multiplet_data = {}
    primary_bc_map = {}
    bc_sets = {}
    for x, y in expanded_data.items():
        jac = y[5]
        if jac >= min_jac:
            multiplet_data[x] = y
            a, b = x

            a_primary = primary_bc_map.setdefault(a, a)
            b_primary = primary_bc_map.setdefault(b, b)
            a_set = bc_sets.setdefault(a_primary, set([a])) # initialize starting sets if needed
            b_set = bc_sets.setdefault(b_primary, set([b]))
            set_info = {a_primary: a_set, b_primary: b_set}

            if a_primary != b_primary:
                remaining_primary = max([a_primary, b_primary], key=barcode_counts.get)
                other_primary = a_primary if remaining_primary == b_primary else b_primary
                for k in set_info[other_primary]:
                    primary_bc_map[k] = remaining_primary

                a_set |= b_set
                bc_sets[remaining_primary] = a_set
                bc_sets.pop(other_primary)

    blacklist = set()
    with open(barcodes_strict, 'w') as f:
        f.write("Barcode1\tBarcode2\tBarcode1Counts\tBarcode2Counts\tCommon\tJaccardIndex\tPrimaryBarcode\n")
        for x, data in multiplet_data.items():
            a, b = x
            pb = primary_bc_map[a]
            if a != pb:
                blacklist.add(a)
            if b != pb:
                blacklist.add(b)
            data[-1] = pb
            f.write("{}\t{}\t{}\t{}\t{}\t{:.4f}\t{}\n".format(*data))

    with open(barcodes_status, 'w') as f:
        f.write("Barcode\tIsMultiplet\tPrimaryBarcode\n")
        for b in barcode_counts.keys():
            if b not in barcodes_considered:
                f.write(f"{b}\tIndeterminate\tNone\n")
            elif b in multiplet_data:
                f.write(f"{b}\tTrue\t{primary_bc_map[b]}\n")
            else:
                f.write(f"{b}\tFalse\tNone\n")

    print_and_log(
        f"Identified {len(multiplet_data)} barcode pairs above Jaccard threshold",
        logout,
        starttime,
    )

    print_and_log(
        f"Identified {len(primary_bc_map)} barcodes belonging to multiplets",
        logout,
        starttime,
    )

    print_and_log(
        f"Identified {len(bc_sets)} unique multiplets",
        logout,
        starttime,
    )

    print_and_log(
        f"After multiplet exclusions, have {len(barcode_counts) - len(blacklist)} total cell barcodes",
        logout,
        starttime,
    )

    logout.close()

if __name__ == '__main__':
    try:
        fragments = sys.argv[1]

        barcodes_strict = sys.argv[2]
        barcodes_expanded = sys.argv[3]
        barcodes_status = sys.argv[4]
        summary = sys.argv[5]
        jac_plot = sys.argv[6]

        main(fragments, barcodes_strict, barcodes_expanded, summary, barcodes_status, jac_plot,min_counts=50, max_beads_per_drop=20, min_cutoff=1e-3)

    except NameError:
        main('/dev/stdin', '/dev/stdout', '/dev/null', '/dev/null', '/dev/null')
