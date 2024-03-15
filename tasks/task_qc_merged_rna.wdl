version 1.0

# TASK
# qc-merged-rna

task qc_merged_rna {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: merged RNA QC task'
    }

    input {
        Array[File] barcode_metadata
        String? prefix
        String? genome_name

        Int? umi_min_cutoff = 1
        Int? gene_min_cutoff = 1
        Int? hist_min_umi = 100
        Int? hist_max_umi = 5000

        # Runtime
        Float? disk_factor = 10.0
        Float? memory_factor = 1.0
        String? docker_image = "us.gcr.io/buenrostro-share-seq/task_qc_rna:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(barcode_metadata, "G")

    # Determining memory size based on the size of the input files.
    Float mem_gb = 4.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String merged_barcode_metadata = "~{prefix}.qc.rna.~{genome_name}.barcode.metadata.tsv"
    String umi_barcode_rank_plot = "~{prefix}.rna.qc.~{genome_name}.umi.barcode.rank.plot.png"
    String gene_barcode_rank_plot = "~{prefix}.rna.qc.~{genome_name}.gene.barcode.rank.plot.png"
    String gene_umi_scatter_plot = "~{prefix}.rna.qc.~{genome_name}.gene.umi.scatter.plot.png"
    String umi_histogram_plot = "~{prefix}.rna.qc~{genome_name}.umi.histogram.png"
    String monitor_log = "monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) 1>&2 &

        # Concatenate barcode metadata files
        echo "------ START: Concatenate barcode metadata files ------" 1>&2
        head -n 1 ~{barcode_metadata[0]} > ~{merged_barcode_metadata}
        for metadata in ~{sep=" " barcode_metadata};
        do
            tail -n +2 $metadata >> ~{merged_barcode_metadata}
        done

        # Make QC plots
        echo "------ START: Generate QC plots ------" 1>&2
        Rscript $(which rna_qc_plots.R) \
            ~{merged_barcode_metadata} \
            ~{umi_min_cutoff} \
            ~{gene_min_cutoff} \
            ~{hist_min_umi} \
            ~{hist_max_umi} \
            ~{umi_barcode_rank_plot} \
            ~{gene_barcode_rank_plot} \
            ~{gene_umi_scatter_plot} \
            ~{umi_histogram_plot}
    >>>

    output {
        File rna_merged_barcode_metadata = "~{merged_barcode_metadata}"
        File? umi_barcode_rank_plot = "~{umi_barcode_rank_plot}"
        File? gene_barcode_rank_plot = "~{gene_barcode_rank_plot}"
        File? gene_umi_scatter_plot = "~{gene_umi_scatter_plot}"
        File? umi_histogram_plot = "~{umi_histogram_plot}"
    }

    runtime {
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        memory: "${mem_gb} GB"
    }
}