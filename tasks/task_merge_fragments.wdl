version 1.0

# TASK
# task-merge-fragments




task atac_merge_fragments {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Merge fragments ATAC'
    }

    input {
        # This task takes in input the aligned bam file and rmeove the low quality reads, the extra chromosomes, marks
        # the duplicats, and convert to a bedpe file.
        Array[File] fragments
        String genome_name
        String? prefix = "sample"
        ## Runtime
        Int? cpus = 2
        Float? disk_factor = 4.0
        Float? memory_factor = 0.08
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_filter_atac"
        String singularity_image = "docker://us.gcr.io/buenrostro-share-seq/share_task_filter_atac"
    }


    # Determine the size of the input
    Float input_file_size_gb = size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 6.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String fragment_file = '${prefix}.atac.final.fragments.${genome_name}.tsv'

    String monitor_log = "atac_merge_fragments_monitor.log"

    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &

        time gzip -dc ~{sep=" " fragments} | sort --parallel ~{cpus} -S 6G -k1,1 -k2,2n -k3,3n | bgzip -c ~{fragments} > ~{fragment_file}.gz

        # ~{prefix}.fragments.tsv.gz.tbi
        tabix --zero-based --preset bed ~{fragment_file}.gz
    >>>

    output {
        File atac_final_fragments = "~{fragment_file}.gz"
        File? atac_final_fragments_index = "~{fragment_file}.gz.tbi"
    }

    runtime {
        cpu: cpus
        memory: "${mem_gb} GB"
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        singularity: "${singularity_image}"
        maxRetries: 1
    }

    parameter_meta {
        fragments: {
                description: 'Fragment file',
                help: 'Array of fragment files to merge',
                example: 'fragments.tsv.gz'
            }
        genome_name: {
                description: 'Reference name.',
                help: 'The name of the reference genome used. This is appended to the output file name.',
                examples: ['GRCh38', 'mm10']
            }
        prefix: {
                description: 'Prefix for output files.',
                help: 'Prefix that will be used to name the output files',
                examples: 'my-experiment'
            }
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus.',
                examples: '4'
            }
        memory_factor: {
                description: 'Multiplication factor to determine memory required for task filter.',
                help: 'This factor will be multiplied to the size of bams to determine required memory of instance (GCP/AWS) or job (HPCs).',
                default: 0.15
            }
        disk_factor: {
                description: 'Multiplication factor to determine disk required for task filter.',
                help: 'This factor will be multiplied to the size of bams to determine required disk of instance (GCP/AWS) or job (HPCs).',
                default: 8.0
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for the filtering step.',
                example: ["us.gcr.io/buenrostro-share-seq/share_task_filter"]
            }
    }


}
