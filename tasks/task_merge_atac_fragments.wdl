version 1.0

# TASK
# merge-atac-fragments

task merge_fragments {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-seq pipeline: merge ATAC fragments task'
    }

    input {
        Array[File] fragments
        String? prefix

        String? docker_image = 'us.gcr.io/buenrostro-share-seq/task_merge_atac_fragments:dev'
        Float? disk_factor = 2
        Int? cpus = 8
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragments, 'G')

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(20 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then 'SSD' else 'LOCAL'

    String output_file = '~{prefix}.fragments.tsv.gz'
    String monitor_log = 'monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &
        
        # decompress fragment files, merge sort, bgzip
        pigz -p ~{cpus} -dc ~{sep=' ' fragments} | cat | sort -k1,1 -k2,2n --parallel=~{cpus} | bgzip -c -@ ~{cpus} > ~{output_file}

        # tabix index
        tabix --zero-based -p bed ~{output_file} 
    >>>

    output {
        File fragments = '${output_file}'
        File fragments_index = '${output_file}.tbi'
        File monitor_log = '${monitor_log}'
    }

    runtime {
        cpu: cpus
        disks: 'local-disk ~{disk_gb} ~{disk_type}'
        docker: '~{docker_image}'
    }

    parameter_meta {
        fragments: {
                description: 'Fragment files',
                help: 'Array of fragment files, one per entity to be merged.',
                example: ['first.fragment.file.tsv.gz', 'second.fragment.file.tsv.gz']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
    }
}
