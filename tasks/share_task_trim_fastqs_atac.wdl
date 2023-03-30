version 1.0

# TASK
# 10x_task_preprocess

task share_trim_fastqs_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: preprocess 10x ATAC data.'
    }

    input {
        # This task takes in input the 3 fastqs coming out from cellranger mkfastqs and preprocess them.
        File fastq_R1 # Pair 1 reads
        File fastq_R2 # Pair 2 reads
        String? prefix

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
    String final_fastq_R1 = "${default="multiome_atac" prefix}.atac.preprocess.cleaned.trimmed.R1.fastq.gz"
    String final_fastq_R2 = "${default="multiome_atac" prefix}.atac.preprocess.cleaned.trimmed.R2.fastq.gz"
    String trimming_log_json = "${default="multiome_atac" prefix}.atac.preprocess.trimming.log.json"
    String trimming_log_html = "${default="multiome_atac" prefix}.atac.preprocess.trimming.log.html"
    String trimming_stats = "${default="multiome_atac" prefix}.atac.preprocess.trimming.adapter.stats.txt"

    String monitor_log = 'monitor_10x_preprocessing.log.txt'

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # read trimming

        fastp -i ~{fastq_R1} -I ~{fastq_R2} -o ~{final_fastq_R1} -O ~{final_fastq_R2} -h ~{trimming_log_html} -j ~{trimming_log_json} -G -Q -L -w ~{cpus} 2> ~{trimming_stats}
    >>>

    output {
        File fastq_R1_trimmed = final_fastq_R1
        File fastq_R2_trimmed = final_fastq_R2
        File tenx_trimming_log_json = trimming_log_json
        File trimming_log_html = trimming_log_html
        File tenx_trimming_stats = trimming_stats
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
