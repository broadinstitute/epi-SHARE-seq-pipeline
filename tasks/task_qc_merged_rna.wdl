version 1.0

# TASK
# qc-merged-rna

task qc_merged_rna {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: QC merged RNA task'
    }

    input {
        File barcode_metadata
        Int? umi_cutoff = 10
        Int? gene_cutoff = 10
        String? genome_name
        String? prefix

        Float? disk_factor = 1.0
        Float? memory_factor = 1.5
        String docker_image = "us.gcr.io/buenrostro-share-seq/task_qc_merged_rna"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(barcode_metadata, "G")

    # Determining memory size based on the size of the input files.
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String umi_barcode_rank_plot = "~{default="merged" prefix}.qc.rna.~{genome_name}.umi.barcode.rank.plot.png"
    String gene_barcode_rank_plot = "~{default="merged" prefix}.qc.rna.~{genome_name}.gene.barcode.rank.plot.png"
    String gene_umi_scatter_plot = "~{default="merged" prefix}.qc.rna.~{genome_name}.gene.umi.scatter.plot.png"
    String monitor_log = "monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Make QC plots
        Rscript $(which rna_qc_plots.R) ~{barcode_metadata} ~{umi_cutoff} ~{gene_cutoff} ~{umi_barcode_rank_plot} ~{gene_barcode_rank_plot} ~{gene_umi_scatter_plot}
    >>>

    output {
        File? umi_barcode_rank_plot = "~{umi_barcode_rank_plot}"
        File? gene_barcode_rank_plot = "~{gene_barcode_rank_plot}"
        File? gene_umi_scatter_plot = "~{gene_umi_scatter_plot}"
    }

    runtime {
        memory : "~{mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ~{disk_gb} ~{disk_type}"
        docker : "${docker_image}"
        maxRetries:1
    }

    parameter_meta {
        barcode_metadata: {
                description: 'Merged RNA barcode metadata file',
                help: 'TSV file containing barcodes and associated numbers of UMIs and genes.',
                example: 'merged_rna_barcode_metadata.tsv'
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
        umi_cutoff: {
                description: 'UMI cutoff',
                help: 'Cutoff for number of UMIs required when making UMI barcode rank plot.',
                example: 10
            }
        gene_cutoff: {
                description: 'Gene cutoff',
                help: 'Cutoff for number of genes required when making gene barcode rank plot.',
                example: 10
            }
    }
}
