version 1.0

# TASK
# SHARE-merge-fragments

task share_merge_fragments {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: merge fragments task'
    }

    input {
        Array[File] fragments
        String? prefix

        String? docker_image = 'us.gcr.io/buenrostro-share-seq/share_task_merge_fragments'
        Float? disk_factor = 1.5
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragments, 'G')

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(5 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then 'SSD' else 'LOCAL'

    String output_file = '${default='aggregated' prefix}.tsv.gz'
    String monitor_log = 'monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        cat ~{sep=' ' fragments} > ~{output_file}
    >>>

    output {
        File fragments = '${output_file}'
        File monitor_log = '${monitor_log}'
    }

    runtime {
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker : "${docker_image}"
    }
}
