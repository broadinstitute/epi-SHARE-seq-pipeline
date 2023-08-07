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
        String chemistry

        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_trim_fastqs_atac:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 16.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # Read trimming outfiles
    String fastq_R1_trimmed = basename(fastq_R1, ".fastq.gz") + "_trimmed.fastq"
    String fastq_R2_trimmed = basename(fastq_R2, ".fastq.gz") + "_trimmed.fastq"
    String trimming_log_json = basename(fastq_R1, "R1.fastq.gz") + ".atac.preprocess.trimming.log.json"
    String trimming_log_html = basename(fastq_R1, "R1.fastq.gz") + ".atac.preprocess.trimming.log.html"
    String trimming_stats = basename(fastq_R1, "R1.fastq.gz") + ".atac.preprocess.trimming.adapter.stats.txt"
    String monitor_log = 'trim_fastqs_atac_monitor.log'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Use trim_fastq script for SHARE ATAC trimming
        if [ '~{chemistry}' == 'shareseq' ]; then
            python3 $(which trim_fastq.py) ~{fastq_R1} ~{fastq_R2} ~{fastq_R1_trimmed} ~{fastq_R2_trimmed} ~{trimming_stats}

        # Use fastp for 10X ATAC trimming
        else
            fastp -i ~{fastq_R1} -I ~{fastq_R2} -o ~{fastq_R1_trimmed} -O ~{fastq_R2_trimmed} -h ~{trimming_log_html} -j ~{trimming_log_json} -G -Q -L -w ~{cpus} 2> ~{trimming_stats}
      
        fi
  
        pigz -p ~{cpus} *.fastq
    >>>

    output {
        File fastq_R1_trimmed = fastq_R1_trimmed + ".gz"
        File fastq_R2_trimmed = fastq_R2_trimmed + ".gz"
        File? tenx_trimming_log_json = trimming_log_json
        File? tenx_trimming_log_html = trimming_log_html
        File trimming_stats = trimming_stats
        File trim_fastqs_atac_monitor = monitor_log
    }

    runtime {
        cpu: cpus
        docker: "${docker_image}"
        disks: "local-disk ${disk_gb} ${disk_type}"
        memory: "${mem_gb} GB"
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
