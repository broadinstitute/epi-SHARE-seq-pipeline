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

chrom_sizes_dict = defaultdict(int)

with open(chrom_sizes, "r") as fh:
    for line in fh:
        chrom_sizes_dict[line.strip().split("\t")[0]] = int(line.strip().split("\t")[1])


data = snap.pp.import_data(
    fragment_file,
    chrom_sizes=chrom_sizes_dict,
    sorted_by_barcode=False,
    min_num_fragments=min_frag_cutoff,
    shift_left=4,
    shift_right=-4
)

snap.metrics.tsse(data, compressed_gtf_file)
snap.pl.tsse(data, min_fragment=min_frag_cutoff, width=800, height=1000, show=False, out_file=tsse_output_file)
snap.metrics.frip(data, {"tss_frac": tss_bed_file, "promoter_frac": promoter_bed_file})

data.obs.to_csv(metrics_output_file, index_label="barcode", sep="\t")
