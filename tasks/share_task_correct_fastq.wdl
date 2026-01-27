version 1.0

# TASK
# SHARE-correct-fastq

task share_correct_fastq {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Correct FASTQs task'
    }

    input {
        File fastq_R1
        File fastq_R2
        File? fastq_barcode
        File whitelist
        String sample_type
        String? pkr
        String? prefix

        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.08
        String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_correct_fastq:v1.0.0"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 16.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String corrected_fastq_R1 = basename(fastq_R1, ".fastq.gz") + "_corrected.fastq"
    String corrected_fastq_R2 = basename(fastq_R2, ".fastq.gz") + "_corrected.fastq"
    String monitor_log = "correct_fastqs_monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Perform barcode error correction on FASTQs
        python3 $(which correct_fastq.py) \
            ~{fastq_R1} \
            ~{fastq_R2} \
            ~{corrected_fastq_R1} \
            ~{corrected_fastq_R2} \
            ~{whitelist} \
            ~{sample_type} \
            ~{prefix} \
            ~{pkr} \
            ~{if defined(fastq_barcode) then "--barcode_fastq " + fastq_barcode else ""}

        pigz -p ~{cpus} *.fastq
    >>>

    output {
        File corrected_fastq_R1 = "~{corrected_fastq_R1}.gz"
        File corrected_fastq_R2 = "~{corrected_fastq_R2}.gz"
        File barcode_qc = "~{prefix}_barcode_qc.txt"
	File monitor_log = "~{monitor_log}"
    }

    runtime {
        cpu : cpus
        memory : "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} ~{disk_type}"
        docker : "~{docker_image}"
    }
}
