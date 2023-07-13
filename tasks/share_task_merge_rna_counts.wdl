version 1.0

# TASK
# SHARE-merge-rna-counts

task share_merge_counts {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: merge RNA counts task'
    }

    input {
        Array[File] tars
        Array[String]? pkrs = ''
        String? gene_naming
        String? prefix

        String? docker_image = 'us.gcr.io/buenrostro-share-seq/share_task_merge_rna_counts'
        Float? disk_factor = 2.0
        Float? memory_factor = 120.0
    }

    # Determine the size of the input
    Float input_file_size_gb = size(tars, 'G')

    # Determining memory size based on the size of the input files.
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then 'SSD' else 'LOCAL'

    String ensembl_option = if '~{gene_naming}'=='ensembl' then '--ensembl' else ''
    String monitor_log = 'monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        if [ '~{sep='' pkrs}' == '' ]; then
            pkr=''
        else
            pkr='--pkr ~{sep=' ' pkrs}'
        fi

        # Create merged h5 matrix
        python3 $(which merge_rna_counts.py) \
            ~{prefix} \
            ~{sep=' ' tars} \
            $pkr \
            ~{ensembl_option} \
    >>>

    output {
        File h5_matrix = '${prefix}.h5'
        File barcode_metadata = '${prefix}_rna_barcode_metadata.tsv'
        File monitor_log = '${monitor_log}'
    }

    runtime {
        memory : "${mem_gb} GB"
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
    }

    parameter_meta {
        tars: {
                description: 'STARsolo output tar.gz files',
                help: 'Array of tar.gz files containing raw matrix, features, and barcodes files from STARsolo, one per entity to be merged.',
                example: ['first.raw.tar.gz', 'second.raw.tar.gz']
            }
        gene_naming: {
                description: 'Gene naming convention',
                help: 'Convention for gene naming in h5 matrix; either "gene_name" (default) or "ensembl".',
                example: ['gene_name', 'ensembl']
            }
    }
}
