import numpy as np
import scipy
import snapatac2 as snap
import sys

from collections import defaultdict

fragment_file = sys.argv[1]
compressed_gtf_file = sys.argv[2]
chrom_sizes = sys.argv[3]
min_frag_cutoff = int(sys.argv[4])
fragments_per_chromosome_file = sys.argv[5]
snap_h5ad = sys.argv[6]
prefix = sys.argv[7]


chrom_sizes_dict = defaultdict(int)

print('Loading chromosomes', file=sys.stderr)
with open(chrom_sizes, "r") as fh:
    for line in fh:
        chrom_sizes_dict[line.strip().split("\t")[0]] = int(line.strip().split("\t")[1])

print('Running import', file=sys.stderr)
data = snap.pp.import_data(
    fragment_file=fragment_file,
    chrom_sizes=chrom_sizes_dict,
    file=snap_h5ad,
    sorted_by_barcode=True,
    min_num_fragments=min_frag_cutoff,
    shift_left=0,
    shift_right=1
)

#data.obs["sample"] = prefix

cell_chr_mat = snap.pp.add_tile_matrix(data, bin_size=1000000000, inplace=False)
cell_chr_mat.write(fragments_per_chromosome_file,  as_dense="X")

data.close()
