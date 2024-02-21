import numpy as np
import scipy
import snapatac2 as snap
import sys

from collections import defaultdict

h5ad = sys.argv[1]
blacklist = sys.argv[2]
chrom_sizes = sys.argv[3]
output = sys.argv[4]

chrom_sizes_dict = defaultdict(int)

print('Loading chromosomes', file=sys.stderr)
with open(chrom_sizes, "r") as fh:
    for line in fh:
        chrom_sizes_dict[line.strip().split("\t")[0]] = int(line.strip().split("\t")[1])


print("Read h5", file=sys.stderr)
counts = snap.read(h5ad, backed=None)
print("Call peaks", file=sys.stderr)
snap.tl.macs3(counts, groupby=["~{prefix}" for i in range(counts.obs["sample"].shape[0])], key_added='peaks', blacklist=blacklist, shift=-75, extsize=150, n_jobs=1)
print("Merge peaks", file=sys.stderr)
peaks = snap.tl.merge_peaks(counts.uns['peaks'], chrom_sizes_dict)
peakcts = snap.pp.make_peak_matrix(counts, use_rep=peaks['Peaks'])
peakcts.write_h5ad(output)

counts.close()
