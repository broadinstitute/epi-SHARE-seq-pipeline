version 1.0

# TASK
# trim_fastqs_atac

task share_trim_fastqs_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: trim ATAC FASTQs.'
    }

    input {
        File fastq_R1 # Pair 1 reads
        File fastq_R2 # Pair 2 reads

        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_trim_fastqs_atac"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # read trimming outfiles
    String fastq_R1_trimmed = basename(fastq_R1, ".fastq.gz") + "trimmed.fastq"
    String fastq_R2_trimmed = basename(fastq_R2, ".fastq.gz") + "trimmed.fastq"

    String monitor_log = 'trim_fastqs_atac_monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # read trimming
        python3 $(which trim_fastq.py) ~{fastq_R1} ~{fastq_R2} ~{fastq_R1_trimmed} ~{fastq_R2_trimmed}
        gzip *.fastq
    >>>

    output {
        File fastq_R1_trimmed = fastq_R1_trimmed + ".gz"
        File fastq_R2_trimmed = fastq_R2_trimmed + ".gz"
        File trim_fastqs_atac_monitor = monitor_log
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        maxRetries: 1
        memory: "${mem_gb} GB"
        memory_retry_multiplier: 2
    }

    parameter_meta {
        fastq_R1: {
                description: 'Pairs 1 fastq',
                help: 'Pairs 1 fastq',
            }
        fastq_R2: {
                description: 'Pairs 2 fastq',
                help: 'Pairs 2 fastq',
            }
    }

}
