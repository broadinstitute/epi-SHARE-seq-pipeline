version 1.0

# TASK
# SHARE-merge-barcode-metadata

task share_merge_barcode_metadata {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: merge barcode metadata task'
    }

    input {
        Array[File] barcode_metadata
        String modality
        String? prefix
        Boolean concat_barcodes = false
        
        String? docker_image = 'us.gcr.io/buenrostro-share-seq/share_task_merge_barcode_metadata'
        Float? disk_factor = 2.0
        Float? memory_factor = 1.0
    }

    # Determine the size of the input
    Float input_file_size_gb = size(barcode_metadata, 'G')

    # Determining memory size based on the size of the input files.
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then 'SSD' else 'LOCAL'

    String output_file = '~{default='aggregated' prefix}.~{modality}.barcode.metadata.tsv'
    String concat_barcodes_option = if concat_barcodes then '--concat' else ''
    String monitor_log = 'monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        python3 $(which merge_barcode_metadata.py) \
            ~{modality} \
            ~{output_file} \
            ~{sep=' ' barcode_metadata} \
            ~{concat_barcodes_option}
    >>>

    output {
        File barcode_metadata = '${output_file}'
        File monitor_log = '${monitor_log}'
    }

    runtime {
        memory : "${mem_gb} GB"
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
    }
}
