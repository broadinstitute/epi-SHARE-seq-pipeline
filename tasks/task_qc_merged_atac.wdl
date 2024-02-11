version 1.0

# TASK
# qc-merged-atac

task qc_merged_atac {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: merged ATAC QC task'
    }

    input {
        Array[File] barcode_metadata
        File fragments
        File fragments_index
        File tss
        Array[String] dataset_names
        String? prefix
        String? genome_name

        Int? fragment_min_cutoff = 10
        Int? hist_max_fragment = 5000

        # Runtime
        Float? disk_factor = 10.0
        Float? memory_factor = 2.0
        String docker_image = "us.gcr.io/buenrostro-share-seq/task_qc_atac:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(barcode_metadata, "G") + size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 4.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String barcode_metadata = '${prefix}.atac.qc.${genome_name}.metadata.tsv'
    String dataset_barcodes = '${prefix}.atac.qc.${genome_name}.dataset.barcodes.tsv'
    String insert_size_hist = '${prefix}.atac.qc.hist.${genome_name}.png'
    String fragment_barcode_rank_plot = '${prefix}.atac.qc.${genome_name}.fragment.barcode.rank.plot.png'
    String fragment_histogram = '${prefix}.atac.qc.${genome_name}.fragment.histogram.png'

    command <<<
        set -e

        bash $(which monitor_script.sh) 1>&2 &

        # Concatenate barcode metadata files
        echo '------ START: Concatenate barcode metadata files ------' 1>&2
        head -n 1 ~{barcode_metadata[0]} > ~{barcode_metadata}
        tail -n +2 ~{sep=' ' barcode_metadata} >> ~{barcode_metadata}

        # Make TSV containing dataset names for each barcode
        echo '------ START: Making dataset barcodes tsv ------' 1>&2
        for i in range(length(~{barcode_metadata}));
        do
            cut -f1 ${barcode_metadata[$i]} | awk -v dataset="${dataset_names[$i]}" -v OFS="\t" 'NR>1{print $0, dataset}' > ~{dataset_barcodes}
        done

        # Insert size plot bulk
        gzip -dc ~{fragments} | awk '{print $3-$2}' > insert_sizes
        echo '------ START: Generate TSS enrichment plot for bulk ------' 1>&2
        time python3 $(which plot_insert_size_hist.py) insert_sizes ~{prefix} ~{insert_size_hist}

        # Barcode rank plot
        echo '------ START: Generate barcode rank plot ------' 1>&2
        time Rscript $(which atac_qc_plots.R) ~{barcode_metadata} ~{fragment_min_cutoff} ~{hist_max_fragment} ~{fragment_barcode_rank_plot} ~{fragment_histogram}
    >>>

    output {
        File atac_barcode_metadata = barcode_metadata
        File atac_dataset_barcodes = dataset_barcodes
        File insert_size_hist = insert_size_hist
        File? fragment_barcode_rank_plot = fragment_barcode_rank_plot
        File? fragment_histogram = fragment_histogram
    }

    runtime {
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        barcode_metadata: {
                description: 'ATAC barcode metadata files',
                help: 'Array of TSV files, each containing barcodes and associated statistics; one per entity to be merged.',
                example: ['first.atac.qc.metadata.tsv', 'second.atac.qc.metadata.tsv']
            }  
        fragments: {
                description: 'Merged fragment file',
                help: 'Merged ATAC fragment file.',
                example: 'merged.fragments.tsv.gz'
            }
        fragments_index: {
                description: 'Merged fragment file index',
                help: 'Index for merged ATAC fragment file.',
                example: 'merged.fragments.tsv.gz.tbi'
            }
        tss: {
                description: 'TSS bed file',
                help: 'List of TSSs in BED format.',
                example: 'refseq.tss.bed'
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
    }    
}
