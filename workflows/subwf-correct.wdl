version 1.0

import "../tasks/share_task_correct_fastq.wdl" as correct_fastq

workflow correct_fastq {
    meta {
        version: 'v0.1'
            author: 'Mei Knudson (mknudson@broadinstitute.org) @ Broad Institute of MIT and Harvard'
            description: 'Broad Institute of MIT and Harvard: Correct FASTQs'
    }

    input {
        Array[File] read1
        Array[File] read2
        File whitelist
        String sample_type
        String? pkr
        String? prefix

        Int? cpus
        Float? disk_factor
        Float? memory_factor
        String? docker_image
    }

    scatter (read_pair in zip(read1, read2)) {
        call correct_fastq.share_correct_fastq as correct {
            input:
                fastq_R1 = read_pair.left,
                fastq_R2 = read_pair.right,
                whitelist = whitelist,
                sample_type = sample_type,
                pkr = pkr,
                prefix = prefix,
                cpus = cpus,
                disk_factor = disk_factor,
                memory_factor = memory_factor,
                docker_image = docker_image
        }
    }

    output {
        Array[File] corrected_fastq_R1 = correct.corrected_fastq_R1
        Array[File] corrected_fastq_R2 = correct.corrected_fastq_R2
        Array[File] corrected_fastq_barcode = correct.corrected_fastq_barcode
    }
}
