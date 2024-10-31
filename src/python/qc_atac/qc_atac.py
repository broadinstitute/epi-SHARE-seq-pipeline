import click
import logging
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import snapatac2 as snap
import sys

from collections import defaultdict
from matplotlib.ticker import FuncFormatter


# Define a custom formatter for the y-axis
def thousands_formatter(x, pos):
    return '%1.0fk' % (x * 1e-3)


# Configure logging
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

@click.command()
@click.option('--fragment_file', type=click.Path(exists=True), help='Path to the fragment file.', required=True)
@click.option('--compressed_gtf_file', type=click.Path(exists=True), help='Path to the compressed GTF file.', required=True)
@click.option('--chrom_sizes', type=click.Path(exists=True), help='Path to the chromosome sizes file.', required=True)
# @click.option('--tss_bed_file', type=click.Path(exists=True), help='Path to the TSS BED file.', required=False)
@click.option('--min_frag_cutoff', default=10, type=int, help='Minimum fragment cutoff. Default is 10.')
@click.option('--prefix', default='snap', type=str, help='Prefix for output files. Default is "snapatac_output".')
def main(fragment_file, compressed_gtf_file, chrom_sizes, min_frag_cutoff, prefix):
    chrom_sizes_dict = defaultdict(int)

    logging.info('Loading chromosomes')
    with open(chrom_sizes, "r") as fh:
        for line in fh:
            chrom_sizes_dict[line.strip().split("\t")[0]] = int(line.strip().split("\t")[1])

    logging.info('Running import')
    data = snap.pp.import_data(
        fragment_file=fragment_file,
        chrom_sizes=chrom_sizes_dict,
        sorted_by_barcode=True,
        min_num_fragments=min_frag_cutoff,
        shift_left=4,
        shift_right=-4,
        file=f"{prefix}_snap.h5ad"
    )

    logging.info('Running metrics')
    snap.pl.frag_size_distr(data, max_recorded_size=1000, show=False)
    logging.info('Running TSS enrichment')
    snap.metrics.tsse(data, compressed_gtf_file)
    
    logging.info('Plot fragment size distribution')
    plt.plot(data.uns['frag_size_distr'])
    plt.xlabel("Fragment size")
    plt.ylabel("Number of fragments")
    plt.gca().yaxis.set_major_formatter(FuncFormatter(thousands_formatter))
    plt.savefig(f"{prefix}_fragment_size_distribution.png")
    plt.close()
    
    logging.info('Plot TSS enrichment')
    plt.plot(data.uns['TSS_profile']/data.uns['TSS_profile'][0])
    plt.xlabel("TSS")
    plt.ylabel("TSS enrichment")
    plt.xticks(ticks=[0, 1000, 2000, 3000, 4000], labels=['-2000', '-1000', '0', '1000', '2000'])
    plt.savefig(f"{prefix}_TSS_enrichment.png")
    plt.close()

    logging.info('Plot fraction of duplicates')
    plt.hist(data.obs['frac_dup'], bins=100)
    plt.xlabel("Fraction of duplicates")
    plt.ylabel("Barcode count")
    plt.savefig(f"{prefix}_fraction_of_duplicates.png")
    plt.close()

    logging.info('Plot fraction of mitochondrial fragments')
    plt.hist(data.obs['frac_mito'], bins=100)
    plt.xlabel("Fraction of mitochondrial fragments")
    plt.ylabel("Barcode count")
    plt.savefig(f"{prefix}_fraction_of_mitochondrial_fragments.png")
    plt.close()

    logging.info('Plot fraction of fragments overlapping TSS')
    with open(f"{prefix}_frac_overlap_TSS.txt","w") as fh:
        fh.write(f"{data.uns['frac_overlap_TSS']:.2f}")
    
    logging.info('Plot library TSS enrichment')
    with open(f"{prefix}_library_TSS.txt","w") as fh:
        fh.write(f"{data.uns['library_tsse']:.2f}")

    logging.info('Plot knee plot')
    fig, ax = plt.subplots(figsize=(10, 7))
    knee = np.sort((np.array(data.obs['n_fragment']))) [::-1]
    cell_set = np.arange(len(knee))
    num_cells = cell_set[::-1][0]

    # Calculate the number of cells with more than 500 fragments
    num_cells_above_500 = np.sum(knee > 500)
    num_cells_above_1000 = np.sum(knee > 1000)

    # plot the knee for the current sublibrary
    ax.loglog(cell_set, knee, linewidth=5, label=f"{num_cells} cells\n{num_cells_above_500} cells > 500 fragments\n{num_cells_above_1000} cells > 1000 fragments")

    ax.axhline(y=500, linewidth=1.5, color="k", linestyle='--')
    ax.axhline(y=1000, linewidth=1.5, color="k", linestyle='--')

    ax.set_ylabel("Number of fragments", fontsize=18)
    ax.set_xlabel("Barcode rank",fontsize=18)
    ax.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')

    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.grid(True, which="both")
    plt.savefig(f"{prefix}_knee_plot.png")
    plt.close()

    logging.info('Plot number of fragments vs TSS enrichment')
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(3, 3, width_ratios=[4, 1, 0.1], height_ratios=[1, 4, 0.1])

    ax_main = plt.subplot(gs[1, 0])
    ax_top = plt.subplot(gs[0, 0], sharex=ax_main)
    ax_right = plt.subplot(gs[1, 1], sharey=ax_main)

    # Main scatter plot
    ax_main.scatter(np.log10(data.obs['n_fragment']), data.obs['tsse'], s=10)
    ax_main.set_xlabel('log10(Number of fragments)')
    ax_main.set_ylabel('TSS enrichment')

    # Add the number of cells in the upper right corner
    num_cells = data.obs['n_fragment'].shape[0]
    ax_main.text(0.95, 0.95, f'Number of cells: {num_cells}', transform=ax_main.transAxes, 
                 verticalalignment='top', horizontalalignment='right', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

    # Top histogram
    ax_top.hist(np.log10(data.obs['n_fragment']), bins=50, color='gray')
    ax_top.axis('off')

    # Right histogram
    ax_right.hist(data.obs['tsse'], bins=100, orientation='horizontal', color='gray')
    ax_right.axis('off')

    plt.tight_layout()
    plt.savefig(f"{prefix}_n_fragment_vs_TSS_enrichment.png", bbox_inches='tight')
    plt.close()

    logging.info('Plot number of fragments vs TSS enrichment (filtered)')
    # Filter the data to keep only barcodes with 500 or more fragments
    filtered_data = data[data.obs['n_fragment'] >= 500]

    # Update the plots to use the filtered data
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(3, 3, width_ratios=[4, 1, 0.1], height_ratios=[1, 4, 0.1])

    ax_main = plt.subplot(gs[1, 0])
    ax_top = plt.subplot(gs[0, 0], sharex=ax_main)
    ax_right = plt.subplot(gs[1, 1], sharey=ax_main)

    # Main scatter plot
    ax_main.scatter(np.log10(filtered_data.obs['n_fragment']), filtered_data.obs['tsse'], s=10)
    ax_main.set_xlabel('log10(Number of fragments)')
    ax_main.set_ylabel('TSS enrichment')

    # Add the number of cells in the upper right corner
    num_cells_filtered = filtered_data.shape[0]
    ax_main.text(0.95, 0.95, f'Number of cells after filtering: {num_cells_filtered}', transform=ax_main.transAxes, 
                verticalalignment='top', horizontalalignment='right', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

    # Top histogram
    ax_top.hist(np.log10(filtered_data.obs['n_fragment']), bins=50, color='gray')
    ax_top.axis('off')

    # Right histogram
    ax_right.hist(filtered_data.obs['tsse'], bins=100, orientation='horizontal', color='gray')
    ax_right.axis('off')

    plt.tight_layout()
    plt.savefig(f"{prefix}_n_fragment_vs_TSS_enrichment_filtered.png", bbox_inches='tight')
    plt.close()

    snap.pp.filter_cells(data, min_counts=500, min_tsse=7, max_counts=100000,)
    snap.pp.add_tile_matrix(data, bin_size=500, counting_strategy="paired-insertion", min_frag_size=20, max_frag_size=1000)

    snap.pp.select_features(data, n_features=200000)

    snap.tl.spectral(data)
    snap.tl.umap(data)
    snap.pp.knn(data)
    snap.tl.leiden(data)
    snap.pl.umap(data, color='leiden', interactive=False, height=500, show=False, out_file=f"{prefix}_umap_leiden.png")

    # Save the content of data.obs to a file
    data.obs.to_csv(f"{prefix}_barcode_metrics.csv")
    logging.info('Done')
    data.close()


if __name__ == '__main__':
    main()
