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

snap_h5ad = sys.argv[10]
prefix = sys.argv[11]


chrom_sizes_dict = defaultdict(int)

print('Loading chromosomes', file=sys.stderr)
with open(chrom_sizes, "r") as fh:
    for line in fh:
        chrom_sizes_dict[line.strip().split("\t")[0]] = int(line.strip().split("\t")[1])

print('Running import', file=sys.stderr)
data = snap.pp.import_data(
    fragment_file=fragment_file,
    chrom_sizes=chrom_sizes_dict,
    filename=f"{prefix}.snapatac.h5ad",
    sorted_by_barcode=True,
    min_num_fragments=100,
    shift_left=0,
    shift_right=1
)
print('TSSe', file=sys.stderr)
snap.metrics.tsse(data, compressed_gtf_file)

print('Plot TSSe', file=sys.stderr)
snap.pl.tsse(data, min_fragment=min_frag_cutoff, width=800, height=1000, show=False, out_file=f"{prefix}.cutoff{min_frag_cutoff}.tsse.png")

print('Plot TSSe strict', file=sys.stderr)
snap.pl.tsse(data, min_fragment=500, width=800, height=1000, show=False, out_file=f"{prefix}.cutoff500.tsse.png")

print('Plot Insert size distribution', file=sys.stderr)
snap.pl.frag_size_distr(data, width=800, height=1000, show=False, out_file=f"{prefix}.cutoff{min_frag_cutoff}.insertsize.distribution.png")

print('Fraction in promoter', file=sys.stderr)
snap.metrics.frip(data, {"tss_frac": tss_bed_file, "promoter_frac": promoter_bed_file})


# Plot frac_mito vs frac_dup
# create usable reads chromap
# create mapped reads
# Frac dup after 500 

print('save metrics', file=sys.stderr)
data.obs.to_csv(f"{prefix}.barcode_stats.tsv", index_label="barcode", sep="\t")

data.obs["sample"] = prefix

snap.pp.filter_cells(data, min_counts=500, min_tsse=0)
cell_chr_mat = snap.pp.add_tile_matrix(data, bin_size=1000000000, inplace=False)
cell_chr_mat.write(f"{prefix}.chromosome_counts_matrix.npz",  as_dense="X")

snap.ex.export_coverage(data,
                        groupby="sample",
                        bin_size=25,
                        blacklist=None,
                        normalization='CPM',
                        min_frag_length=None,
                        max_frag_length=2000,
                        out_dir='./', prefix=prefix,
                        suffix='.bw',
                        tempdir='./',
                        n_jobs=8)

data.close()
