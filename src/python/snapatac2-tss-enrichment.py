import numpy as np
import scipy
import snapatac2 as snap
import sys

from collections import defaultdict

fragment_file = sys.argv[1]
compressed_gtf_file = sys.argv[2]
chrom_sizes = sys.argv[3]
tss_bed_file = sys.argv[4]
promoter_bed_file = sys.argv[5]
min_frag_cutoff = int(sys.argv[6])
metrics_output_file = sys.argv[7]
tsse_output_file = sys.argv[8]
fragments_per_chromosome_file = sys.argv[9]


chrom_sizes_dict = defaultdict(int)

print('Loading chromosomes', file=sys.stderr)
with open(chrom_sizes, "r") as fh:
    for line in fh:
        chrom_sizes_dict[line.strip().split("\t")[0]] = int(line.strip().split("\t")[1])

print('Running import', file=sys.stderr)
data = snap.pp.import_data(
    fragment_file=fragment_file,
    chrom_sizes=chrom_sizes_dict,
    sorted_by_barcode=True,
    min_num_fragments=min_frag_cutoff,
    shift_left=0,
    shift_right=1
)
print('TSSe', file=sys.stderr)
snap.metrics.tsse(data, compressed_gtf_file)
print('Plot TSSe', file=sys.stderr)
snap.pl.tsse(data, min_fragment=min_frag_cutoff, width=800, height=1000, show=False, out_file=tsse_output_file)
print('Fraction in promoter', file=sys.stderr)
snap.metrics.frip(data, {"tss_frac": tss_bed_file, "promoter_frac": promoter_bed_file})
print('save metrics', file=sys.stderr)
data.obs.to_csv(metrics_output_file, index_label="barcode", sep="\t")

snap.pp.filter_cells(data, min_counts=np.quantile(data.obs["n_fragment"], q=[.9])[0], min_tsse=0)
cell_chr_mat = snap.pp.add_tile_matrix(data, bin_size=1000000000, inplace=False)
scipy.sparse.save_npz(file=fragments_per_chromosome_file, cell_chr_mat.X)
