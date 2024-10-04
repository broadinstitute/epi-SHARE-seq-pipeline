version 1.0

import "../tasks/share_task_correct_fastq.wdl" as correct_fastq

workflow correct_fastq {
    meta {
        version: 'v0.1'
            author: 'Mei Knudson (mknudson@broadinstitute.org) @ Broad Institute of MIT and Harvard'
            description: 'Broad Institute of MIT and Harvard: Correct FASTQs'
    }

    input {
        Array[File] atac_fastq_R1
        Array[File] atac_fastq_R2
        Array[File] rna_fastq_R1
        Array[File] rna_fastq_R2
        File whitelist
        String? pkr
        String? prefix
        
        # Runtime attributes
        Int? cpus
        Float? disk_factor
        Float? memory_factor
        String? docker_image
    }

    if ( length(atac_fastq_R1) > 0 ) {
        scatter (read_pair in zip(atac_fastq_R1, atac_fastq_R2)) {
            call correct_fastq.share_correct_fastq as correct_atac {
                input:
                    fastq_R1 = read_pair.left,
                    fastq_R2 = read_pair.right,
                    whitelist = whitelist,
                    sample_type = "ATAC",
                    pkr = pkr,
                    prefix = prefix,
                    cpus = cpus,
                    disk_factor = disk_factor,
                    memory_factor = memory_factor,
                    docker_image = docker_image
            }
        }
    }

    if ( length(rna_fastq_R1) > 0 ) {
        scatter (read_pair in zip(rna_fastq_R1, rna_fastq_R2)) {
            call correct_fastq.share_correct_fastq as correct_rna {
                input:
                    fastq_R1 = read_pair.left,
                    fastq_R2 = read_pair.right,
                    whitelist = whitelist,
                    sample_type = "RNA",
                    pkr = pkr,
                    prefix = prefix,
                    cpus = cpus,
                    disk_factor = disk_factor,
                    memory_factor = memory_factor,
                    docker_image = docker_image
            }
        }
    }

    output {
        Array[File]? atac_corrected_fastq_R1 = correct_atac.corrected_fastq_R1
        Array[File]? atac_corrected_fastq_R2 = correct_atac.corrected_fastq_R2
        Array[File]? atac_corrected_fastq_barcode = correct_atac.corrected_fastq_barcode

        Array[File]? rna_corrected_fastq_R1 = correct_rna.corrected_fastq_R1
        Array[File]? rna_corrected_fastq_R2 = correct_rna.corrected_fastq_R2
    }
}
